from pycule import MCuleWrapper

import os


class MCuleAPI(object):
    """
    Creates an MCuleAPI object for calling Mcule API
    """

    def __init__(self):
        """
        MCule API constructor
        """
        self.mculewrapper = MCuleWrapper(authorisation_token=os.environ["MCULE_API_KEY"])

    def getMCuleInfo(self, smiles: str):
        """
        Get compound info from Mcule

        Args:
            smiles (str): SMILES string of compound to search MCule DB for
        Returns:
            mculeid: MCule ID for compound
            mculeurl: MCule url for compound
            None: If no inof is found, returns None
        """
        try:
            response_dict = self.mculewrapper.singlequerysearch(query=smiles)
            results = response_dict["response"]["results"]
            if results:
                mculeid = results[0]["mcule_id"]
                mculeurl = results[0]["url"]
                return mculeid, mculeurl
            else:
                return None
        except Exception as e:
            print(e)

    def getMCulePrice(self, mculeid: str, amount: int = 10):
        """
        Get compound pricing info from Mcule for 1, 5 and 10 mg amounts

        Args:
            mculeid (str): MCule ID for compound
            amount (int): 1, 5 or 10 mg for gettig price estimate for
        Returns:
            price: MCule price in USD for comppound and amount set in iput argument
        """
        try:
            amount_dict = {1: 0, 5: 1, 10: 2}
            response_dict = self.mculewrapper.compoundprices(mcule_id=mculeid)
            bestprices = response_dict["response"]["best_prices"]
            amountpriceinfo = bestprices[amount_dict[amount]]
            if amountpriceinfo:
                price = amountpriceinfo["price"]
                deliverytime = amountpriceinfo["delivery_time_working_days"]
                return price, deliverytime
            else:
                return None
        except Exception as e:
            print(e)
