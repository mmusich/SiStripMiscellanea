Run the file inspector.
1) prepare the library:
```
root -b 
[0].L printSizes.C++g
[1].q
```
2) run the file macro to extract the json file:
```
python singleFileInspector.py INPUT_FILE OUTPUT_JSON_FILE
```

notice that the `INPUT_FILE` name should be in the format `calibTree_<run>.root`.
Then run the script to compare the two outputfiles:

```
python pretty_print.py output_new.json output_old.json
```