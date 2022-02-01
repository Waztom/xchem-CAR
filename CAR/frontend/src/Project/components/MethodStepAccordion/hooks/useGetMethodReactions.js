import axios from 'axios';
import { useEffect, useState } from 'react';

export const useGetMethodReactions = (methods) => {
  const [methodReactions, setMethodReactions] = useState([]);

  useEffect(() => {
    const fetch = async () => {
      const requests = methods.map(async (method) => {
        const response = await axios.get(`/api/reactions?search=${method.id}`);
        return {
          method,
          reactions: response.data,
        };
      });
      const responses = await Promise.allSettled(requests);

      setMethodReactions(
        responses
          .filter((response) => response.status === 'fulfilled')
          .map((response) => response.value)
      );

      responses
        .filter((response) => response.status === 'rejected')
        .forEach((response) => console.error(response.reason));
    };

    fetch();
  }, [methods]);

  return methodReactions;
};
