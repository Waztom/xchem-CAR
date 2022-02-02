export const useCategorizeMethodsByNoSteps = (methodsWithTarget) => {
  const categorizedMethodsWithTarget = {};
  methodsWithTarget.forEach((methodWithTarget) => {
    const noSteps = methodWithTarget.method.nosteps;

    let methodsForStep = categorizedMethodsWithTarget[noSteps];
    if (!methodsForStep) {
      methodsForStep = [];
      categorizedMethodsWithTarget[noSteps] = methodsForStep;
    }

    methodsForStep.push(methodWithTarget);
  });

  return categorizedMethodsWithTarget;
};
