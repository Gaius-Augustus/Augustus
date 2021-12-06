# AUGUSTUS Command Line Parameters
All possible command line parameters of AUGUSTUS are listed in the JSON file `aug_cmdln_parameters.json`. In addition to the actual name, the file also contains the following additional information about the respective parameter.

## Keys for JSON configuration file
| Key | Type | Description |
| --- | ---- | ----------- |
| name |string | The name of the parameter. |
| type | string | The type of the parameter value. Possible types: `string`, `int`, `float`, `bool`, `list<string>`. |
| possible_values | [string] | A list of possible values for the parameter. |
| usage | string | Examples of use of the parameter. |
| default_value | string | A default value for the parameter. |
| development | bool | If this key is set to true, the parameter is not yet intended for productive use or information is still missing. |
| exclude_apps | string | Programs for which this parameter is not valid. Possible values are: `augustus`, `etraining` and `pygustus`.|
| description | string | The description of the parameter. Should descriptions be created for a list of possible values, key-value pairs can also be specified here.  As an example for this the parameter `genemodel` can be considered.
