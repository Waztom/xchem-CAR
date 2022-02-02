import { makeStyles } from '@material-ui/core';

const useStyles = makeStyles((theme) => ({
  icon: {
    width: 24,
    height: 24,
  },
}));

export const IconComponent = ({ Component }) => {
  const classes = useStyles();

  return <Component className={classes.icon} />;
};
