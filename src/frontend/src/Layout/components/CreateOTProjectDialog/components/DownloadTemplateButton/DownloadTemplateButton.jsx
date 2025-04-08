import React from 'react';
import { Button } from '@mui/material';
import { Download } from '@mui/icons-material';

export const DownloadTemplateButton = () => {
  const handleDownload = () => {
    const csvContent = `labware_type,well_index,smiles,volume,concentration,solvent
fluidx_24_vials_2500ul,0,CC(=O)O,100,0.1,DMSO
fluidx_96_vials_1000ul,1,CCO,200,0.2,MeOH`;

    const blob = new Blob([csvContent], { type: 'text/csv' });
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'starting_materials_template.csv';
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    window.URL.revokeObjectURL(url);
  };

  return (
    <Button
      startIcon={<Download />}
      onClick={handleDownload}
      size="small"
    >
      Download Template
    </Button>
  );
};