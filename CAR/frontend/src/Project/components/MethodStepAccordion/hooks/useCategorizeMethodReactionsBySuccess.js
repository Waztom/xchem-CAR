export const useCategorizeMethodReactionsBySuccess = (methodReactions) => {
  const categorizedMethodReactions = {};

  methodReactions.forEach((methodReaction) => {
    const successfulReactions = methodReaction.reactions.filter(
      (reaction) => reaction.successrate >= 0.5
    ).length;

    let methodReactionsForSuccess =
      categorizedMethodReactions[successfulReactions];
    if (!methodReactionsForSuccess) {
      methodReactionsForSuccess = [];
      categorizedMethodReactions[successfulReactions] =
        methodReactionsForSuccess;
    }

    methodReactionsForSuccess.push(methodReaction);
  });

  return categorizedMethodReactions;
};
