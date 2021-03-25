# document of Synthetic Photometry(synphot)

#download of data to utilize the pre-defined standard star, extinction laws, and bandpasses
from synphot.utils import download_data
file_list = download_data('./synphot_data')
