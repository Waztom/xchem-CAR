import { MenuItem } from '@material-ui/core';
import React from 'react';
import { useGetOtProjects } from '../../../../../common/hooks/useGetOtProjects';
import { useCurrentProjectStore } from '../../../../../common/stores/currentProjectStore';

export const OTProjectMenuContents = ({ onSelected }) => {
  const currentProject = useCurrentProjectStore.useCurrentProject();

  const { data: otProjects } = useGetOtProjects({ project_id: currentProject.id });

  return !!otProjects.length ? (
    otProjects.map(otProject => {
      return (
        <MenuItem key={otProject.id} onClick={() => onSelected(otProject)}>
          {otProject.name}
        </MenuItem>
      );
    })
  ) : (
    <MenuItem disabled>There are no OT protocols</MenuItem>
  );
};
