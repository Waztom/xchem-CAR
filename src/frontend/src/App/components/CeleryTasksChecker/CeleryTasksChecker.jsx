import React from 'react';
import { useCeleryTasksStore } from '../../../common/stores/celeryTasksStore';
import { CeleryTask } from './components/CeleryTask/CeleryTask';

export const CeleryTasksChecker = () => {
  const tasks = useCeleryTasksStore.useTasks();

  return (
    <>
      {Object.values(tasks).map(task => (
        <CeleryTask key={task.id} task={task} />
      ))}
    </>
  );
};
