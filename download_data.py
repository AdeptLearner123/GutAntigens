import os
import requests

url = 'https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/raw/HMP2/MGX/2018-05-04/HSMA33OT.tar'

def main():
    os.mkdir("genomes", exist_ok=True)
    print ("Sending request")
    response = requests.get(url)
    print ("Received response")
    if response.status_code == 200:
        print ("Writing")
        with open("genomes/test", 'wb') as file:
            file.write(response.content)
    print('File downloaded successfully')


if __name__ == "__main__":
    main()