The scripts in this folder are used to convert DB2 files to db2bin format.
The `db2_to_binary_recursive.py` has been tested to convert DB2.tgz files in the ZINC22 database to db2bin.tar.zst format (it uses zstandard compression)

The performance gain is about ~17x faster than text parsing and the file size is about 20% smaller than the original DB2.tgz files.

To use the binary format in DOCK3, you need to specify in the INDOCK:

```
INDOCK: use_binary_db2: yes
```