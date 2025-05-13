import os
import requests
import tarfile

url = 'https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/raw/HMP2/MGX/2018-05-04/HSMA33OT.tar'

def main():
    """
    os.makedirs("genomes", exist_ok=True)
    print ("Sending request")
    response = requests.get(url)
    print ("Received response")
    if response.status_code == 200:
        print ("Writing")
        with open("genomes/test", 'wb') as file:
            file.write(response.content)
    print('File downloaded successfully')
    """
    tar = tarfile.open("genomes/test", "r:tar")
    tar.extractall("destination_folder")
    tar.close()


if __name__ == "__main__":
    main()