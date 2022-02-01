import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableRow,
} from '@material-ui/core';

export const ReactionTable = ({ noSteps, methodData }) => {
  return (
    <Table>
      <TableHead>
        <TableRow>
          <TableCell>Compound</TableCell>
          {new Array(noSteps).fill(0).map((_, index) => {
            return <TableCell key={index}>Reaction {index + 1}</TableCell>;
          })}
        </TableRow>
      </TableHead>
      <TableBody>
        {methodData.map(({ method, reactions, target }) => {
          return (
            <TableRow key={method.id}>
              <TableCell>
                <img src={target.image} height={50} />
              </TableCell>
              {reactions.map((reaction) => {
                return (
                  <TableCell key={reaction.id}>
                    <img src={reaction.reactionimage} height={50} />
                  </TableCell>
                );
              })}
            </TableRow>
          );
        })}
      </TableBody>
    </Table>
  );
};
