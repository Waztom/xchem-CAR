import { useEffect, useState } from 'react';
import axios from 'axios';

export const useGetCategorizedMethodsBySteps = (projectId) => {
  const [methods, setMethods] = useState([]);

  useEffect(() => {
    if (!projectId) {
      return;
    }

    async function fetchData() {
      try {
        const response = await axios.get(`/api/methods/`);
        setMethods(response.data);
      } catch (err) {
        console.error(err);
      }
    }
    fetchData();
  }, [projectId]);

  const sortMethods = () => {
    const sortedMethods = {};
    methods.forEach((method) => {
      const noSteps = method.nosteps;

      let methodsForStep = sortedMethods[noSteps];
      if (!methodsForStep) {
        methodsForStep = [];
        sortedMethods[noSteps] = methodsForStep;
      }

      methodsForStep.push(method);
    });

    return sortedMethods;
  };

  return sortMethods();
};
