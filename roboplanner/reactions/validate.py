from rdkit import Chem

def add_warning(field, warning_string, validate_dict):
    validate_dict['field'].append(field)
    validate_dict['warning_string'].append(warning_string)

    return validate_dict

def checkColumnNames(columns, validate_dict):
    column_names = ['Targets', 'Ammount_required (mg)']
    if not all(columns == column_names):
        validate_dict = add_warning(field = 'name_columns',
                                    warning_string = 'Column names should be set to: {}'.format(column_names),
                                    validate_dict=validate_dict)
    return validate_dict


def checkNumberColumns(columns, validate_dict):
    no_columns = len(columns)
    
    if no_columns > 2:
        validate_dict = add_warning(field = 'number_columns',
                                    warning_string = "Found {} column names. Set and name columns to 'Targets' only".format(no_columns),
                                    validate_dict=validate_dict)
        
    if no_columns == 2:
        validate_dict = checkColumnNames(columns, validate_dict)
    
    return validate_dict


def checkSMILES(target_smiles, index, validate_dict):
    
    mol = Chem.MolFromSmiles(target_smiles)
    
    if mol is None:
        validate_dict = add_warning(field = 'check_smiles',
                                    warning_string = "Input target smiles: '{}' at index {} is not a valid smiles".format(target_smiles, index),
                                    validate_dict=validate_dict)

    return validate_dict

def checkIsNumber(amount, index, validate_dict):
    if type(amount) != float:
         validate_dict = add_warning(target_name = target_smiles,
                            field = 'check_number',
                            warning_string = "Input target smiles: '{}' at index {} is not a valid number".format(target_smiles, index),
                            validate_dict=validate_dict)
    
    return validate_dict
