from pycule import MCuleWrapper
import os
import inspect

from ..utils import canonSmiles

import logging

logger = logging.getLogger(__name__)


class MCuleAPI(object):
    """
    Creates an MCuleAPI object for calling Mcule API
    """

    def __init__(self):
        """
        MCule API constructor
        """
        self.mculewrapper = MCuleWrapper(
            authorisation_token=os.environ["MCULE_API_KEY"]
        )

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

        except Exception as e:
            logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
            try:
                response_dict = self.mculewrapper.similaritysearch(
                    query=smiles, limit=1, threshold=0.7
                )
                results = response_dict["response"]["results"]
                if results:
                    smiles_test = canonSmiles(results[0]["smiles"])
                    if smiles_test == smiles:
                        mculeid = results[0]["mcule_id"]
                        mculeurl = results[0]["url"]
                        return mculeid, mculeurl
                    else:
                        return None
            except Exception as e:
                logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
                print(e)
                return None

    def getMCulePrice(self, mculeid: str, amount: float):
        """
        Get compound pricing info from Mcule for 1, 5 and 10 mg amounts

        Args:
            mculeid (str): MCule ID for compound
            amount (floatt): Amount required
        Returns:
            price: MCule price in USD for comppound and amount set in iput argument
        """
        if amount < 1:
            amount = 1
        try:
            response_dict = self.mculewrapper.compoundpricesamount(
                mcule_id=mculeid, amount=amount
            )
            price_info = response_dict["response"]["best_prices"][0]
            if price_info:
                price = price_info["price"]
                deliverytime = price_info["delivery_time_working_days"]
                return price, deliverytime

        except Exception as e:
            logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
            print(e)
            print(response_dict)
            return None

    def getTotalQuote(
        self,
        mculeids: list,
        amount: float = 1,
        delivery_country: str = "GB",
    ):
        """
        Get quote from MCule for list of mcule ids

        Args:
            mculeids (list): List of MCule IDs
            amount (float): Amount per compound for quote
            delivery_country (str): ISO 3166-1 alpha-2 code of the delivery country. Default GB
            target_volume (float): Total volume in ml requested. Default None
            target_cc (float): Target concentration in mM. Default None
        Returns:
            quote (dict): MCule quote info as a dictionary
        """
        quote_info = {}

        try:
            response_dict = self.mculewrapper.quoterequest(
                mcule_ids=mculeids,
                delivery_country=delivery_country,
                amount=amount,
            )
            quote_id = response_dict["response"]["id"]
            quote_state_response = self.mculewrapper.quoterequeststatus(
                quote_id=quote_id
            )
            quote_state = quote_state_response["response"]["state"]

            while quote_state != 30:
                if quote_state == 40:
                    return None
                else:
                    quote_state_response = self.mculewrapper.quoterequeststatus(
                        quote_id=quote_id
                    )
                    quote_state = quote_state_response["response"]["state"]

            quote_info["quoteid"] = quote_state_response["response"]["group"]["quotes"][
                0
            ]["id"]
            quote_info["quoteurl"] = quote_state_response["response"]["group"][
                "quotes"
            ][0]["site_url"]
            quote_info["quotecost"] = quote_state_response["response"]["group"][
                "quotes"
            ][0]["total_cost"]
            quote_info["quotevaliduntil"] = quote_state_response["response"]["group"][
                "quotes"
            ][0]["valid_until"]

            return quote_info

        except Exception as e:
            logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
            print(e)
