# OpenbindDatabaseTemplate

This contains the format of the tables used to train the PanDDA 2 event scoring models and the code required to construct them from XChem data.

For performance reasons tables are stored as zarr arrays rather than in sql (some table items are very large and it was considered better to keep to one format for data rather than two). The datatypes in `src/xchem_dataset/tables` give the tabular descriptions of the data, and the foriegn key structure should be fairly appararent from variable names. 

The code has been layed out flat, as many of the functions in here are frequently useful for working with XChem data!