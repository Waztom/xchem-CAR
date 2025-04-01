import React from 'react';
import { Typography, Table, TableBody, TableRow, TableCell } from '@mui/material';
import { styled } from '@mui/material/styles';
import {
  setErrorsForSmilesValidationErrorsDialog,
  setSmilesValidationErrorsDialogOpen,
  useSmilesValidationErrorsDialogStore
} from '../../stores/smilesValidationErrorsDialogStore';
import { CloseDialog } from '../../../common/components/CloseDialog';
import { DialogSection } from '../../../common/components/DialogSection';
import { DialogSectionHeading } from '../../../common/components/DialogSectionHeading';

const StyledTableRow = styled(TableRow)(({ theme }) => ({
  '&:nth-of-type(odd)': {
    backgroundColor: theme.palette.action.hover
  }
}));

export const SmilesValidationErrorsDialog = () => {
  const { dialogOpen, errors } = useSmilesValidationErrorsDialogStore();

  return (
    <CloseDialog
      id="smiles-validation-errors-dialog"
      open={dialogOpen}
      title="Smiles validation errors"
      content={
        <DialogSection>
          <DialogSectionHeading>Errors</DialogSectionHeading>
          <Typography>This is a list of errors which were encountered during smiles validation:</Typography>
          <Table>
            <TableBody>
              {errors?.warning_string.map((error, index) => (
                <StyledTableRow key={index}>
                  <TableCell>
                    {index + 1}. {error}
                  </TableCell>
                </StyledTableRow>
              ))}
            </TableBody>
          </Table>
        </DialogSection>
      }
      onClose={() => {
        setSmilesValidationErrorsDialogOpen(false);
      }}
      TransitionProps={{
        onExited: () => setErrorsForSmilesValidationErrorsDialog(null)
      }}
    />
  );
};
