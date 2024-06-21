import os
import numpy as np
from keras.preprocessing import image
from keras.applications.resnet50 import ResNet50, preprocess_input
from sklearn.cluster import KMeans
from keras.models import Model
import pandas as pd
import pickle

# Constants
BATCH_SIZE = 50
TARGET_SIZE = (224, 224)


# Load and preprocess images in batches
def load_images_in_batches(folder, batch_size=BATCH_SIZE):
    batch_images = []
    batch_filenames = []
    for filename in os.listdir(folder):
        try:
            img_path = os.path.join(folder, filename)
            img = image.load_img(img_path, target_size=TARGET_SIZE)
            img_array = image.img_to_array(img)
            img_array = preprocess_input(img_array)
            batch_images.append(img_array)
            batch_filenames.append(filename)
            if len(batch_images) == batch_size:
                yield np.array(batch_images), batch_filenames
                batch_images, batch_filenames = [], []  # Reset for next batch
        except Exception as e:
            print(f"Error loading image: {filename}. Error: {e}")
    if batch_images:  # Yield any remaining images
        yield np.array(batch_images), batch_filenames


# Load ResNet50 model
base_model = ResNet50(weights="imagenet")
model = Model(inputs=base_model.input, outputs=base_model.get_layer("avg_pool").output)


# Extract features
def extract_features(images):
    return model.predict(images, batch_size=10)


# Cluster features
def cluster_features(features, n_clusters=5):
    kmeans = KMeans(n_clusters=n_clusters, random_state=0)
    kmeans.fit(features)
    return kmeans.labels_


# Main process
folder_path = "/work/rwth1209/projects/hackathon_ghent/cell_images/"  # Update this to your folder path
all_filenames = []
all_labels = []
all_features = []

for images, filenames in load_images_in_batches(folder_path):
    features = extract_features(images)
    labels = cluster_features(features)
    all_filenames.extend(filenames)
    all_labels.extend(labels)
    all_features.extend(features)

# Save results
results = pd.DataFrame(
    {"filename": all_filenames, "label": all_labels, "features": list(all_features)}
)
# results.to_csv("resnet50_n5_results.csv", index=False)

with open("resnet50_n5_results.pkl", "wb") as f:
    pickle.dump(results, f)
