import axios from 'axios';
import { useEffect, useState } from 'react';

export const useGetTargets = (projectId) => {
  const [targets, setTargets] = useState([]);

  useEffect(() => {
    if (!projectId) {
      return;
    }

    async function fetchData() {
      try {
        const response = await axios.get(`/api/targets?search=${projectId}`);
        setTargets(response.data);
      } catch (err) {
        console.error(err);
      }
    }
    fetchData();
  }, [projectId]);

  return targets;
};
