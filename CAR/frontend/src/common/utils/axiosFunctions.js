import axios from 'axios';

export const axiosGet = async (...params) => {
  const response = await axios.get(...params);
  return response.data;
};

export const axiosPatch = async (...params) => {
  const response = await axios.patch(...params);
  return response.data;
};
