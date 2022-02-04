export const useCategorizeMethodsDataBySuccessRate = (methodReactions) => {
  const categorizedMethodReactions = {};

  /**
   * The way the categorization functions is that for each reaction a success is either 0 or 1 depending on the
   * value of successrate. These 0 and 1 values are then concatenated into a string, which is then easily sorted
   * by the consuming component. They key itself can be split to get the original values.
   */
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
