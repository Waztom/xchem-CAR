import React from 'react';
import { SubmitDialog } from '../../../common/components/SubmitDialog';
import { DialogSection } from '../../../common/components/DialogSection';
import { DialogSectionHeading } from '../../../common/components/DialogSectionHeading';
import { Typography } from '@material-ui/core';
import { Formik, Form } from 'formik';
import * as yup from 'yup';
import { useUploadProject } from './hooks/useUploadProject';
import { FormTextField } from '../../../common/components/FormTextField';
import { FormRadioGroup } from '../../../common/components/FormRadioGroup';
import { FormFilePicker } from '../../../common/components/FormFilePicker';
import { useValidateSmiles } from './hooks/useValidateSmiles';

const validationOptions = [
  { value: '0', label: 'Validate' },
  { value: '1', label: 'Upload' }
];

const apiOptions = [
  { value: '0', label: 'Postera' },
  { value: '1', label: 'Custom chemistry' },
  { value: '2', label: 'Combi custom chemistry' }
];

export const UploadProjectDialog = ({ open, onClose }) => {
  const { mutate: uploadProject } = useUploadProject();
  const { mutate: validateSmiles } = useValidateSmiles();

  return (
    <Formik
      initialValues={{
        project_name: '',
        submitter_name: '',
        submitter_organisation: '',
        protein_target: '',
        csv_file: null,
        validate_choice: '',
        API_choice: ''
      }}
      validationSchema={yup.object().shape({
        project_name: yup
          .string()
          .min(1, 'The name is too short')
          .max(100, 'The name is too long')
          .matches(/^[a-zA-Z0-9\-_ ]+$/, 'The name can only include a-z, A-Z, 0-9, -, _ or space')
          .required('Required'),
        submitter_name: yup
          .string()
          .min(1, 'The name is too short')
          .max(100, 'The name is too long')
          .matches(/^[a-zA-Z0-9\-_ ]+$/, 'The name can only include a-z, A-Z, 0-9, -, _ or space')
          .required('Required'),
        submitter_organisation: yup
          .string()
          .min(1, 'The name is too short')
          .max(100, 'The name is too long')
          .matches(/^[a-zA-Z0-9\-_ ]+$/, 'The name can only include a-z, A-Z, 0-9, -, _ or space')
          .required('Required'),
        protein_target: yup
          .string()
          .min(1, 'The name is too short')
          .max(100, 'The name is too long')
          .matches(/^[a-zA-Z0-9\-_ ]+$/, 'The name can only include a-z, A-Z, 0-9, -, _ or space')
          .required('Required'),
        csv_file: yup.mixed().required('Required'),
        validate_choice: yup
          .string()
          .oneOf(validationOptions.map(({ value }) => value))
          .required('Select one of the choices'),
        API_choice: yup
          .string()
          .oneOf(apiOptions.map(({ value }) => value))
          .required('Select one of the choices')
      })}
      onSubmit={data => {
        if (data.validate_choice === '0') {
          validateSmiles({ data });
        } else {
          uploadProject({ data });
        }
        onClose();
      }}
    >
      {({ isSubmitting, resetForm }) => (
        <SubmitDialog
          id="upload-project-dialog"
          open={open}
          title="Upload new project"
          content={
            <Form id="upload-project-form">
              <DialogSection>
                <DialogSectionHeading>Project information</DialogSectionHeading>
                <Typography>Please provide following information:</Typography>

                <FormTextField
                  name="project_name"
                  label="Project name"
                  placeholder="Project name (can include a-z, A-Z, 0-9, -, _ or space)"
                />

                <FormTextField
                  name="submitter_name"
                  label="Your name"
                  placeholder="Your name (can include a-z, A-Z, 0-9, -, _ or space)"
                />

                <FormTextField
                  name="submitter_organisation"
                  label="Your organisation"
                  placeholder="Your organisation (can include a-z, A-Z, 0-9, -, _ or space)"
                />

                <FormTextField
                  name="protein_target"
                  label="Protein target name"
                  placeholder="Protein target name (can include a-z, A-Z, 0-9, -, _ or space)"
                />

                <FormFilePicker
                  name="csv_file"
                  label="Smiles file"
                  description={
                    <Typography>
                      The current specification version is <strong>ver_1.2</strong>.
                    </Typography>
                  }
                  id="upload-project-smiles-file"
                  accept="text/csv"
                />

                <FormRadioGroup name="validate_choice" label="Validate choice" options={validationOptions} />

                <FormRadioGroup name="API_choice" label="API choice" options={apiOptions} />
              </DialogSection>
            </Form>
          }
          onClose={() => {
            onClose();
          }}
          submitDisabled={isSubmitting}
          SubmitButtonProps={{
            type: 'submit',
            form: 'upload-project-form'
          }}
          TransitionProps={{
            onExited: () => {
              resetForm();
            }
          }}
        />
      )}
    </Formik>
  );
};
