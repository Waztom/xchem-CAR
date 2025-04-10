import React from 'react';
import { Button } from '@mui/material';
import { Download } from '@mui/icons-material';

export const DownloadTemplateButton = () => {
  const handleDownload = () => {
    const csvContent = `plate-ID,labware-type,well-index,SMILES,amount-uL,concentration,solvent
plate1,fluidx_24_vials_2500ul,0,CC(=O)O,100,0.1,DMSO
plate1,fluidx_24_vials_2500ul,1,CC(=O)OCCC,100,0.1,DMSO
plate2,fluidx_96_vials_1000ul,A1,CCO,200,0.2,MeOH
plate2,fluidx_96_vials_1000ul,B1,CCCCO,200,0.2,MeOH`;

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