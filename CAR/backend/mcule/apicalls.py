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

    def getTotalQuote(
        self,
        mculeids: list,
        delivery_country: str = "GB",
        target_volume: float = None,
        target_cc: float = None,
    ):
        """
        Get quote from MCule for list of mcule ids

        Args:
            mculeids (list): List of MCule IDs
            delivery_country (str): ISO 3166-1 alpha-2 code of the delivery country. Default GB
            target_volume (float): Total volume in ml requested. Default None
            target_cc (float): Target concentration in mM. Default None
        Returns:
            quote: MCule quote
        """
        try:
            response_dict = self.mculewrapper.quoterequest(
                mcule_ids=mculeids,
                delivery_country=delivery_country,
                target_volume=target_volume,
                target_cc=target_cc,
            )
            quote_id = response_dict["response"]["id"]
            quote_state_response = self.mculewrapper.quoterequeststatus(quote_id=quote_id)
            quote_state = quote_state_response["response"]["state"]

            while quote_state != 30:
                if quote_state == 40:
                    return None
                else:
                    quote_state_response = self.mculewrapper.quoterequeststatus(quote_id=quote_id)
                    quote_state = quote_state_response["response"]["state"]

            quote_url = quote_state_response["response"]["site_url"]
            quote_cost = quote_state_response["response"]["group"]["quotes"][0]["total_cost"]

            return quote_url, quote_cost

        except Exception as e:
            print(e)
