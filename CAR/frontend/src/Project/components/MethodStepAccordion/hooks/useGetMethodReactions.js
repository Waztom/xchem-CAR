import axios from 'axios';
import { useEffect, useState } from 'react';

export const useGetMethodReactions = (methodsWithTarget) => {
  const [methodsData, setMethodsData] = useState([]);

  useEffect(() => {
    const fetch = async () => {
      const requests = methodsWithTarget.map(async (methodWithTarget) => {
        const response = await axios.get(
          `/api/reactions?search=${methodWithTarget.method.id}`
        );
        return {
          ...methodWithTarget,
          reactions: response.data,
        };
      });
      const responses = await Promise.allSettled(requests);

      setMethodsData(
        responses
          .filter((response) => response.status === 'fulfilled')
          .map((response) => response.value)
      );

      responses
        .filter((response) => response.status === 'rejected')
        .forEach((response) => console.error(response.reason));
    };

    fetch();
  }, [methodsWithTarget]);

  return methodsData;
};
