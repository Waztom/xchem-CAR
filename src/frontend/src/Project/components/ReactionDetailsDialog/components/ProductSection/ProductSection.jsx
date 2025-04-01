import React from 'react';
import { Box, Tooltip, Typography } from '@mui/material';
import { styled } from '@mui/material/styles';
import { DialogSection } from '../../../../../common/components/DialogSection';
import { DialogSectionHeading } from '../../../../../common/components/DialogSectionHeading';
import { SuspenseWithBoundary } from '../../../../../common/components/SuspenseWithBoundary';
import PubChem from '../../../../../assets/pubchem.svg';
import PubChemSafety from '../../../../../assets/pubchem-safety.svg';

const SectionsGrid = styled('div')(({ theme }) => ({
  display: 'grid',
  gridTemplateColumns: 'repeat(3, 1fr)'
}));

const ProductSectionContent = ({ product }) => {
  const productpubcheminfo = product.productpubcheminfo || {};

  return (
    <DialogSection>
      <DialogSectionHeading>Product</DialogSectionHeading>
      <SectionsGrid>
        <Box>
          <Typography>
            Smiles: <strong>{product.smiles}</strong>
          </Typography>
          {productpubcheminfo.cas && (
            <Typography>
              CAS number: <strong>{productpubcheminfo.cas}</strong>
            </Typography>
          )}
        </Box>

        {!!productpubcheminfo.summaryurl && (
          <Box>
            <Tooltip title="PubChem compound summary">
              <a href={productpubcheminfo.summaryurl} target="_blank" rel="noreferrer">
                <img src={PubChem} height={50} alt="PubChem" />
              </a>
            </Tooltip>
          </Box>
        )}
        {!!productpubcheminfo.lcssurl && (
          <Box>
            <Tooltip title="Laboratory chemical safety summary">
              <a href={productpubcheminfo.lcssurl} target="_blank" rel="noreferrer">
                <img src={PubChemSafety} height={50} alt="PubChem Safety" />
              </a>
            </Tooltip>
          </Box>
        )}
      </SectionsGrid>
    </DialogSection>
  );
};

export const ProductSection = (props) => (
  <SuspenseWithBoundary>
    <ProductSectionContent {...props} />
  </SuspenseWithBoundary>
);

ProductSection.displayName = 'ProductSection';
