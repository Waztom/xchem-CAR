from rdkit import Chem

def add_warning(target_name, field, warning_string, validate_dict):
    validate_dict['target_name'].append(target_name)
    validate_dict['field'].append(field)
    validate_dict['warning_string'].append(warning_string)

    return validate_dict

def checkColumnName(column_name, validate_dict):
    if column_name != 'Targets':
        validate_dict = add_warning(target_name = 'Column name error',
                                    field = 'column_name',
                                    warning_string = "Column name set to {} and should be 'Targets'".format(column_name),
                                    validate_dict=validate_dict)
    return validate_dict


def checkNumberColumns(columns, validate_dict):
    no_columns = len(columns)
    
    if no_columns > 1:
        validate_dict = add_warning(target_name = 'Column name error',
                                    field = 'column_name',
                                    warning_string = "Found {} column names. Set and name columns to 'Targets' only".format(no_columns),
                                    validate_dict=validate_dict)
        
    if no_columns == 1:
        validate_dict = checkColumnName(columns[0], validate_dict)
    
    return validate_dict


def checkSMILES(target_smiles, index, validate_dict):
    
    mol = Chem.MolFromSmiles(target_smiles)
    
    if mol is None:
        validate_dict = add_warning(target_name = target_smiles,
                            field = 'smiles_check',
                            warning_string = "Input target smiles: '{}' at index {} is not a valid smiles".format(target_smiles, index),
                            validate_dict=validate_dict)

    return validate_dict