import spatialdata_plot
from dotenv import load_dotenv
from spatialdata import read_zarr
from upath import UPath

load_dotenv()

s3_endpoint_url = "https://s3.embl.de"

bucket_name = "spatialdata/spatialdata-sandbox"
zarr_store_path = "xenium_rep1_io.zarr"

access_key = None  # os.environ.get("ACCESS_KEY")
secret_key = None  # os.environ.get("SECRET_KEY")

# Complete S3 path
s3_path = f"s3://{bucket_name}/{zarr_store_path}"

# Create a UPath object for the .zarr folder with authentication
# Change anon to False for private buckets, and provide access_key and secrey_key
path = UPath(s3_path, anon=True, key=access_key, secret=secret_key, client_kwargs={"endpoint_url": s3_endpoint_url})

if path.exists():
    print(f"The path {s3_path} exists.")
    for item in path.iterdir():
        print(item)
else:
    print(f"The path {s3_path} does not exist.")

sdata = read_zarr(store=path, selection=["images"])

arr = sdata.images["morphology_focus"]["scale0"]["image"].data[0, 1000:1010, 1000:1010].compute()

assert ~(arr == 0).all(), "Data does not seem to be loaded correctly"

# will save figure as test.png in folder
sdata.pp.get_elements(["morphology_focus"]).pl.render_images().pl.show( save = "test" )

# To write locally, do:

# sdata.write( path )