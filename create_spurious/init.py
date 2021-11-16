import os

if not os.path.exists("fake_files"):
    os.makedirs("fake_files")
if not os.path.exists("blastp_search"):
    os.makedirs("blastp_search")
    os.makedirs("blastp_search/blast")
    os.makedirs("blastp_search/queries")
if not os.path.exists("logs"):
    os.makedirs("logs")
