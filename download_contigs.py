import os
import requests
import gzip
import shutil
from tqdm import tqdm
from constants import *

def main():
    print("Fetching URLs...")
    with open(HMP2_CONTIGS_URLS_FILE) as f:
        urls = f.read().splitlines()

    print("Downloading...")
    os.makedirs(CONTIGS_DOWNLOADS_DIR, exist_ok=True)
    for url in tqdm(urls):
        file_name = url.split("/")[-1]
        file_path = os.path.join(CONTIGS_DOWNLOADS_DIR, file_name)
        
        if not os.path.exists(file_path):
            response = requests.get(url)
            if response.status_code == 200:
                with open(file_path, "wb") as f:
                    f.write(response.content)
            else:
                print("Failed to download file at url: " + url + " " + str(response.status_code) + " " + response.text)
                break

    print("Decompressing...")
    files = os.listdir(CONTIGS_DOWNLOADS_DIR)
    os.makedirs(CONTIGS_DIR, exist_ok=True)
    for file_name in tqdm(files):
        file_path = os.path.join(CONTIGS_DOWNLOADS_DIR, file_name)
        file_out_path = os.path.join(CONTIGS_DIR, file_name.replace(".gz", ""))
        if not os.path.exists(file_out_path):
            with gzip.open(file_path, "rb") as f_in:
                with open(file_out_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)


if __name__ == "__main__":
    main()