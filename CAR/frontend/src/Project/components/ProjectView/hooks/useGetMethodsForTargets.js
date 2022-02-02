import axios from 'axios';
import { useEffect, useState } from 'react';

export const useGetMethodsForTargets = (targets) => {
  const [methodsWithTarget, setMethodsWithTarget] = useState([]);

  useEffect(() => {
    const fetch = async () => {
      const requests = targets.map(async (target) => {
        const response = await axios.get(`/api/methods?search=${target.id}`);

        return response.data
          .filter((method) => method.target_id === target.id)
          .map((method) => ({
            target,
            method,
          }));
      });
      const responses = await Promise.allSettled(requests);

      setMethodsWithTarget(
        responses
          .filter((response) => response.status === 'fulfilled')
          .map((response) => response.value)
          .flat()
      );

      responses
        .filter((response) => response.status === 'rejected')
        .forEach((response) => console.error(response.reason));
    };

    fetch();
  }, [targets]);

  return methodsWithTarget;
};
