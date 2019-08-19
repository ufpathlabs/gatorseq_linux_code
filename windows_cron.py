import os
import time

t_end = time.time() + 60 * 60 * 12

while time.time() < t_end:

    print("running QCI_upload.py")
    os.system("C:/ProgramData/Anaconda3/python.exe c:/Users/s.majety/Desktop/gatorseq_linux_code/QCI_windows_download.py")
    time.sleep(10)

    print("running QCI_download.py")
    os.system("C:/ProgramData/Anaconda3/python.exe c:/Users/s.majety/Desktop/gatorseq_linux_code/QCI_windows_upload.py")
    time.sleep(10)

    print("running EPIC_upload.py")
    os.system("C:/ProgramData/Anaconda3/python.exe c:/Users/s.majety/Desktop/gatorseq_linux_code/EPIC_windows_upload.py")
    time.sleep(60)