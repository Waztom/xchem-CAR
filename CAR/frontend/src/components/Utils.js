import axios from "axios";

export function isInt(n) {
  return Number(n) === n && n % 1 === 0;
}

export function isFloat(n) {
  return Number(n) === n && n % 1 !== 0;
}

export async function patchChange(actiontype, id, variablename, value) {
  try {
    const response = await axios.patch(`api/IBM${actiontype}actions/${id}/`, {
      [variablename]: value,
    });
    console.log(response.data);
  } catch (error) {
    console.log(error);
  }
}
