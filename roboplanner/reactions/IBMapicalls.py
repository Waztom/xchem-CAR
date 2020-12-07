from rxn4chemistry import RXN4ChemistryWrapper

# Setup IBM RxN API
# api_key=os.environ['IBM_API_KEY'] 
# rxn4chemistry_wrapper = RXN4ChemistryWrapper(api_key=api_key)
# rxn4chemistry_wrapper.create_project('Test actions')

def getIBMRetroSyn(smiles):
        """
        Use the IBM API to get some possible retrosynthesis routes
        """
        # Create dummy dictionary to create while loop to catch when status is a SUCCESS
        results = {}
        results['status'] = None
        reactions = None

        while reactions is None:  
            try:
                time.sleep(30)
                response = rxn4chemistry_wrapper.predict_automatic_retrosynthesis(product=smiles)
                while results['status'] != 'SUCCESS': 
                    time.sleep(30)
                    results = rxn4chemistry_wrapper.get_predict_automatic_retrosynthesis_results(response['prediction_id'])
                    reactions = results
            except Exception as e:
                print(e)
        return reactions