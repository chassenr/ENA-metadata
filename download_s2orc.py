import requests
import os
from urllib.parse import urlparse, unquote
from tqdm import tqdm
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

API_KEY = "dbQgWRAUPU5FQ82XZlbZNaYTs364Nz8x9i0vrESc"
BASE_URL = "https://api.semanticscholar.org/datasets/v1/release/latest/dataset/s2orc"
DOWNLOAD_DIR = "s2orc_data"
MAX_WORKERS = 4 

os.makedirs(DOWNLOAD_DIR, exist_ok=True)
headers = {"x-api-key": API_KEY}

# Fetch file list
response = requests.get(BASE_URL, headers=headers)
response.raise_for_status()
file_urls = response.json().get("files", [])

def download_file(url, retries=5):
    parsed_url = urlparse(url)
    file_name = unquote(os.path.basename(parsed_url.path))
    file_path = os.path.join(DOWNLOAD_DIR, file_name)
    temp_file_path = file_path + ".part"

    # Check if partially downloaded
    resume_header = {}
    if os.path.exists(temp_file_path):
        resume_header = {"Range": f"bytes={os.path.getsize(temp_file_path)}-"}

    for attempt in range(retries):
        try:
            with requests.get(url, headers={**headers, **resume_header}, stream=True, timeout=30) as r:
                r.raise_for_status()
                mode = "ab" if resume_header else "wb"
                total_size = int(r.headers.get('Content-Length', 0))
                # Adjust total size for resume
                if 'Content-Range' in r.headers:
                    total_size = int(r.headers['Content-Range'].split('/')[-1])
                with open(temp_file_path, mode) as f, tqdm(
                    total=total_size, initial=os.path.getsize(temp_file_path) if os.path.exists(temp_file_path) else 0,
                    unit='B', unit_scale=True, desc=file_name
                ) as pbar:
                    for chunk in r.iter_content(chunk_size=1024*1024):
                        if chunk:
                            f.write(chunk)
                            pbar.update(len(chunk))
            os.rename(temp_file_path, file_path)
            return f"Downloaded: {file_name}"
        except (requests.ConnectionError, requests.Timeout) as e:
            wait = 5 * (attempt + 1)
            print(f"{file_name} failed, retrying in {wait}s ({attempt+1}/{retries})")
            time.sleep(wait)
    return f"Failed: {file_name}"

# Parallel download
with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
    future_to_url = {executor.submit(download_file, url): url for url in file_urls}
    for future in as_completed(future_to_url):
        print(future.result())
