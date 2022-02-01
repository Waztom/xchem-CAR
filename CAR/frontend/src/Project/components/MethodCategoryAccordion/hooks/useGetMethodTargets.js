import axios from 'axios';
import { useEffect, useState } from 'react';

export const useGetMethodTargets = (methodReactions) => {
  const [methodData, setMethodData] = useState([]);

  useEffect(() => {
    const fetch = async () => {
      const requests = methodReactions.map(async (methodReaction) => {
        const response = await axios.get(
          `/api/targets/${methodReaction.method.target_id}`
        );
        return {
          ...methodReaction,
          target: response.data,
        };
      });
      const responses = await Promise.allSettled(requests);

      setMethodData(
        responses
          .filter((response) => response.status === 'fulfilled')
          .map((response) => response.value)
      );

      responses
        .filter((response) => response.status === 'rejected')
        .forEach((response) => console.error(response.reason));
    };

    fetch();
  }, [methodReactions]);

  return methodData;
};
