import PIL.Image
from uni import get_encoder
import torch
import PIL
import os
import pickle
import warnings
from torch.utils.data import DataLoader, Dataset

warnings.filterwarnings("ignore")

IMG_FOLDER = "/work/rwth1209/projects/hackathon_ghent/cell_images"

device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"Using device: {device}")

model, transform = get_encoder(enc_name="uni", device=device)
model.to(device)
model.eval()


class ImageDataset(Dataset):
    def __init__(self, img_folder, transform=None):
        self.img_folder = img_folder
        self.file_names = [
            f
            for f in os.listdir(img_folder)
            if os.path.isfile(os.path.join(img_folder, f))
        ]
        self.transform = transform

    def __len__(self):
        return len(self.file_names)

    def __getitem__(self, idx):
        img_path = os.path.join(self.img_folder, self.file_names[idx])
        try:
            with PIL.Image.open(img_path) as img:
                if self.transform:
                    img = self.transform(img)
        except (IOError, OSError) as e:
            print(f"Error opening image {img_path}: {e}")
            return None, self.file_names[idx]  # Return None for the image
        return img, self.file_names[idx]


# Custom collate function to filter out None values and separate images and filenames
def custom_collate(batch):
    batch = [(img, fname) for img, fname in batch if img is not None]
    if len(batch) == 0:
        return (
            torch.Tensor(),
            [],
        )  # Return empty tensor and list if all images were None
    imgs, fnames = zip(*batch)
    imgs = torch.stack(imgs)  # Stack images into a tensor
    return imgs, list(fnames)


dataset = ImageDataset(IMG_FOLDER, transform=transform)
data_loader = DataLoader(
    dataset, batch_size=50, num_workers=4, pin_memory=True, collate_fn=custom_collate
)

filenames = []
features = []

for images, batch_filenames in data_loader:
    if images.nelement() == 0:  # Check if the images tensor is empty
        continue
    print(f"Processing batch with images: {batch_filenames}")
    images = images.to(device)
    with torch.inference_mode():
        batch_features = model(images)  # Extracted features
    filenames.extend(batch_filenames)
    features.extend(
        batch_features.cpu().detach().numpy()
    )  # Move tensors to CPU and convert to numpy

    # Free GPU memory
    del images
    del batch_features
    torch.cuda.empty_cache()

# Make a dictionary
res = dict(zip(filenames, features))

with open("uni_features.pkl", "wb") as f:
    pickle.dump(res, f)
