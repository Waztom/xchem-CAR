import { useContext } from 'react';
import { BatchContext } from '../context/BatchContext';

export const useBatchContext = () => {
  return useContext(BatchContext);
};
