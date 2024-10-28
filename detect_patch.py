import os
import cv2
import numpy as np
import pandas as pd
import shapely.geometry as geom
import shutil
from tqdm import tqdm
import warnings

# Constants
PATH_ELLIPSES = r"C:\Users\Agathe\Mon Drive\ellipses\\"
PATH_RECO = r"C:\Users\Agathe\Desktop\recolorize_output2\\"
PATH_ANTPOST = r"C:\Users\Agathe\Mon Drive\remove_tail_translate_DV\without_tail_antpost\\"
PATH_SAVE_PATCHES = r"C:\Users\Agathe\Mon Drive\ellipses\000patches\\"
PATH_LAND_F = r"C:\Users\Agathe\Mon Drive\Ariane\Morphométrie\Data brutes\Forewings\Brutes"
PATH_LAND = r"C:\Users\Agathe\Mon Drive\Ariane\Morphométrie\Data brutes\Hindwings\Brutes"
PATCHES_DATA_FILE = r"C:\Users\Agathe\Mon Drive\Données\Morphology\patches_df.csv"

# Clean up directories
def clean_directory(directory):
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print(f'Failed to delete {file_path}. Reason: {e}')

# Calculate the perimeter of an ellipse
def calculate_perimeter(a, b):
    return 2 * np.pi * np.sqrt((a**2 + b**2) / 2)

# Find the image with the largest non-white area (background)
def find_background(images):
    max_pixels = 0
    bg_index = 0
    for i, im_path in enumerate(images):
        im = cv2.imread(im_path, cv2.IMREAD_GRAYSCALE)
        _, bw = cv2.threshold(im, 127, 255, cv2.THRESH_BINARY)
        npix = np.sum(bw == 255)
        if npix > max_pixels:
            max_pixels = npix
            bg_index = i
    return bg_index

# Check if a contour is irregular based on its number of vertices
def is_irregular_shape(contour, vertex_threshold=10):
    epsilon = 0.02 * cv2.arcLength(contour, True)
    approx_polygon = cv2.approxPolyDP(contour, epsilon, True)
    return len(approx_polygon) > vertex_threshold

# Filter out irregular shapes from a binary mask
def filter_irregular_shapes(binary_mask):
    contours, _ = cv2.findContours(binary_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    filtered_contours = [c for c in contours if not is_irregular_shape(c)]
    filtered_mask = np.zeros_like(binary_mask)
    cv2.drawContours(filtered_mask, filtered_contours, -1, 255, cv2.FILLED)
    return filtered_mask

# Setup blob detector
def setup_blob_detector():
    params = cv2.SimpleBlobDetector_Params()
    params.filterByColor = False
    params.filterByArea = False
    params.filterByCircularity = True
    params.minCircularity = 0.5
    params.maxCircularity = 1
    params.filterByConvexity = True
    params.minConvexity = 0.5
    params.maxConvexity = 1
    params.filterByInertia = True
    params.minInertiaRatio = 0.2
    params.maxInertiaRatio = 1
    ver = (cv2.__version__).split('.')
    if int(ver[0]) < 3:
        return cv2.SimpleBlobDetector(params)
    return cv2.SimpleBlobDetector_create(params)

# Load and process landmarks
def load_landmarks(name_im, df, image_ori):
    df = df.dropna()
    df["n"] = np.cumsum(df["X"].notnull() & df["Y"].notnull())
    
    if df["X"].dtype != "float":
        df["X"] = df["X"].str.replace(",", ".").astype(float)
    if df["Y"].dtype != "float":
        df["Y"] = df["Y"].str.replace(",", ".").astype(float)

    df["X"], df["Y"] = df["Y"], df["X"]

    width_scale = image_ori.shape[0] / 2848
    height_scale = image_ori.shape[1] / 4288

    df['X'] *= width_scale
    df['Y'] = image_ori.shape[1] - (df['Y'] * height_scale)

    if name_im[9] == "P":
        if len(df["X"]) == 164:
            CU2 = geom.Point(df[df['n'] == 8]["Y"], df[df['n'] == 8]["X"])
            df = df[df['n'] > 19]
        else:
            print("Landmark error")
            return None, None
    else:
        CU2 = None
        df = df[df['n'] > 18]

    return df, CU2

# Extract patches and save results
def extract_and_save_patches(CU2, A_wing, name_im, image_ori, contours_poly, detector_patch, list_layers_impath, image_ori_big):
    patches_data = []
    
    polygon = geom.Polygon(contours_poly)
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(3,3))
    kernel_blur_size = 5

    for path in list_layers_impath:
        layer = cv2.imread(path)
        _, bw = cv2.threshold(layer, 127, 255, cv2.THRESH_BINARY)
        bw = cv2.morphologyEx(bw, cv2.MORPH_OPEN, kernel)
        bw = cv2.morphologyEx(bw, cv2.MORPH_CLOSE, kernel)
        blurImg = cv2.blur(bw, (kernel_blur_size, kernel_blur_size))
        _, bw2 = cv2.threshold(blurImg, 100, 255, cv2.THRESH_BINARY)
        closing = cv2.morphologyEx(bw2, cv2.MORPH_OPEN, kernel)
        closing = cv2.morphologyEx(closing, cv2.MORPH_CLOSE, kernel)
        cv2.drawContours(closing, [contours_poly], -1, 0, 5)
        
        contours_all, hierarchy = cv2.findContours(closing[:,:,0],cv2.RETR_LIST, cv2.CHAIN_APPROX_NONE)
        
        keypoints = detector_patch.detect(closing)
        minEllipse = []

        for keypoint in keypoints:
            center = (int(keypoint.pt[0]), int(keypoint.pt[1]))
            size = int(keypoint.size / 2)
            min_distance = float('inf')
            closest_contour = None

            for contour in contours_all:
                distance = cv2.pointPolygonTest(contour, center, True)
                if distance >= 0 and distance < min_distance:
                    min_distance = distance
                    closest_contour = contour

            if closest_contour is not None and len(closest_contour) >= 5:
                ellipse = cv2.fitEllipse(closest_contour)
                minEllipse.append(ellipse)

        drawing = image_ori.copy()
        count_patch = 1
        for ellipse in minEllipse:
            A = np.pi / 4 * ellipse[1][1] * ellipse[1][0]
            if A < A_wing / 10 and A > A_wing / 500:
                cv2.ellipse(drawing, ellipse, (255, 0, 255), 3)
                center, axes, angle = ellipse
                point = geom.Point(center)
                distance = point.distance(polygon)
                if distance == 0:
                    distance = polygon.exterior.distance(point)
                distanceCU2 = point.distance(CU2)
                row_data = {'id': name_im, 'x': center[0], 'y': center[1], 'minor_axis' : axes[0], 'major_axis' : axes[1], "distance_contours" : distance, "distance_CU2" : distanceCU2, "A_wing" : A_wing}
                patches_data.append(row_data)

                scaling = image_ori_big.shape[0]/image_ori.shape[0]
                mask = np.zeros_like(image_ori_big)
                mask = cv2.ellipse(mask, (((center[0]*scaling), center[1]*scaling), (int((axes[0]*scaling) * 1.5), int((axes[1]*scaling) * 1.5)), angle), color=(255, 255, 255), thickness=-1)
                result = np.bitwise_and(image_ori_big, mask)
                non_black_pixels = (result != [0, 0, 0])
                non_black_indices = np.any(non_black_pixels, axis=-1)
                y, x = np.where(non_black_indices)
                x_min, x_max = min(x), max(x)
                y_min, y_max = min(y), max(y)
                result = result[y_min:y_max+1, x_min:x_max+1]
                result[np.all(result == [0, 0, 0], axis=-1)] = 255

                cv2.imwrite(f"{PATH_SAVE_PATCHES}{name_im}_{count_patch}.jpg", result)
                count_patch += 1

        cv2.imwrite(f"{PATH_ELLIPSES}{name_im}.jpg", drawing)

    return patches_data

# Main processing function
def main():
    clean_directory(PATH_ELLIPSES)
    detector_patch = setup_blob_detector()

    list_name_im = os.listdir(PATH_RECO)
    list_name_im.remove("cm")
    # List of filenames in PATH_ANTPOST without their ".jpg" extension
    antpost_files = {y.replace(".jpg", "") for y in os.listdir(PATH_ANTPOST)}
    
    # Filter list_name_im based on the presence in the antpost_files set
    list_name_im = [x for x in list_name_im if x in antpost_files]
    listdir_land_F = os.listdir(PATH_LAND_F)
    listdir_land = os.listdir(PATH_LAND)

    names_land_F = [name for name in listdir_land_F if name.endswith(".tps")]
    names_land = [name for name in listdir_land if name.endswith(".tps")]

    patches_data = []

    # Suppressing warnings within the tqdm loop
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # Ignore all warnings
        for name_im in tqdm(list_name_im):
            print(name_im)
            path_image = os.path.join(PATH_ANTPOST, name_im+".jpg")
            image_ori_big = cv2.imread(path_image)
            scale_percent = 17.55618 # percent of original size
            width_ori = int(image_ori_big.shape[1] * scale_percent / 100)
            height_ori = int(image_ori_big.shape[0] * scale_percent / 100)
            dim = (width_ori, height_ori)
            # resize image
            image_ori = cv2.resize(image_ori_big.copy(), dim, interpolation = cv2.INTER_AREA)
    
            list_layers_impath = [os.path.join(PATH_RECO, name_im, filename) for filename in os.listdir(os.path.join(PATH_RECO, name_im))]
            index = find_background(list_layers_impath)
            list_layers_impath.pop(index)
    
            name_lands = [name for name in (names_land_F if name_im[9] == "A" else names_land) if name_im.split("_")[0][:-1] in name.replace("-","")]
            if len(name_lands) != 1:
                print(f"Landmark error for {name_im}")
                continue
            
            name_land = name_lands[0]
            path_land = os.path.join(PATH_LAND_F if name_im[9] == "A" else PATH_LAND, name_land)
            
            df = pd.read_csv(path_land, sep=" ", names=["X", "Y"])
            df, CU2 = load_landmarks(name_im, df, image_ori)
    
    
    
            # Convert the image to grayscale
            im_gray = cv2.cvtColor(image_ori, cv2.COLOR_BGR2GRAY)
            # Apply a binary threshold to create a binary image
            _, thresh_im_cm = cv2.threshold(im_gray, 254, 255, cv2.THRESH_BINARY_INV)
            # Find contours in the binary image
            contours, _ = cv2.findContours(thresh_im_cm, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            contours = max(contours, key = cv2.contourArea)
            
            w,h,x,y = cv2.boundingRect(contours)
            A_wing =cv2.countNonZero(thresh_im_cm)
            
            contours_poly = np.squeeze(contours)
    
            patch_data = extract_and_save_patches(CU2, A_wing, name_im, image_ori, contours_poly, detector_patch, list_layers_impath, image_ori_big)
            patches_data.extend(patch_data)

    df_patches = pd.DataFrame(patches_data)
    df_patches.to_csv(PATCHES_DATA_FILE, index=False, sep=";")

if __name__ == "__main__":
    main()
