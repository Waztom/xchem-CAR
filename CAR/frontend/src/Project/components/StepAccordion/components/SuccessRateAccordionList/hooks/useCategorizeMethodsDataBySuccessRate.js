export const useCategorizeMethodsDataBySuccessRate = (methodReactions) => {
  const categorizedMethodReactions = {};

  methodReactions.forEach((methodReaction) => {
    const key = methodReaction.reactions.reduce(
      (key, reaction) => key + (reaction.successrate >= 0.5 ? '1' : '0'),
      ''
    );

    let methodReactionsForKey = categorizedMethodReactions[key];
    if (!methodReactionsForKey) {
      methodReactionsForKey = [];
      categorizedMethodReactions[key] = methodReactionsForKey;
    }

    methodReactionsForKey.push(methodReaction);
  });

  return categorizedMethodReactions;
};
