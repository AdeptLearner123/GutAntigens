import os
import requests
import tarfile
import gzip
import shutil
from tqdm import tqdm
from constants import *

def download_urls(urls_file, download_dir):
    print(f"Downloading {urls_file}...")
    with open(urls_file) as f:
        urls = f.read().splitlines()
    
    os.makedirs(download_dir, exist_ok=True)
    for url in tqdm(urls):
        file_name = url.split("/")[-1]
        file_path = os.path.join(download_dir, file_name)
        
        if not os.path.exists(file_path):
            response = requests.get(url)
            if response.status_code == 200:
                with open(file_path, "wb") as f:
                    f.write(response.content)
            else:
                print("Failed to download file at url: " + url + " " + str(response.status_code) + " " + response.text)
                break


def decompress_gene_calls():
    print("Decompressing gene calls...")
    files = os.listdir(GENE_CALLS_DOWNLOADS_DIR)
    os.makedirs(GENE_CALLS_DIR, exist_ok=True)
    for file_name in tqdm(files):
        file_path = os.path.join(GENE_CALLS_DOWNLOADS_DIR, file_name)
        tar = tarfile.open(file_path, "r:*")
        tar.extractall(GENE_CALLS_DIR)
        tar.close()


def decompress_contigs():
    print("Decompressing contigs...")
    files = os.listdir(CONTIGS_DOWNLOADS_DIR)
    os.makedirs(CONTIGS_DIR, exist_ok=True)
    for file_name in tqdm(files):
        file_path = os.path.join(CONTIGS_DOWNLOADS_DIR, file_name)
        file_out_path = os.path.join(CONTIGS_DIR, file_name.replace(".gz", ""))
        if not os.path.exists(file_out_path):
            with gzip.open(file_path, "rb") as f_in:
                with open(file_out_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)


def main():
    download_urls(HMP2_GENE_CALLS_URLS_FILE, GENE_CALLS_DOWNLOADS_DIR)
    #download_urls(HMP2_PILOT_GENE_CALLS_URLS_FILE, GENE_CALLS_DOWNLOADS_DIR)
    # ^ Pilot genomic data only contains peptides
    decompress_gene_calls()
    download_urls(HMP2_CONTIGS_URLS_FILE, CONTIGS_DOWNLOADS_DIR)
    #download_urls(HMP2_PILOT_CONTIGS_URLS_FILE, GENE_CALLS_DOWNLOADS_DIR)
    decompress_contigs()


if __name__ == "__main__":
    main()