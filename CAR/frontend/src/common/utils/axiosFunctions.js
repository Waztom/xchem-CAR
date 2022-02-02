import axios from 'axios';

export const axiosGet = async (...params) => {
  const response = await axios.get(...params);
  return response.data;
};
