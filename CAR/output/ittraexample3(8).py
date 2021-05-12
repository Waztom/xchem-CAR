from opentrons import protocol_api
# ittraexample3 for "Example Author" produced by XChem Car (https://car.xchem.diamond.ac.uk)
# metadata
 metadata = {
	'protocolName': 'ittraexample3',
	'author': 'Example Author',
	'description': '('example description',)',
	'apiLevel': '2.9'}

def run(protocol: protocol_api.ProtocolContext):

	# labware
	OrderPlate_1 = protocol.load_labwear('plate_12', '1')
	ReactionPlate_2 = protocol.load_labwear('plate_12', '2')
	tips_96_10_3 = protocol.load_labwear('opentrons_96_filtertiprack_10ul', '3')

	# pipettes
	left_pipette = protocol.load_instrument('p10_single', 'left', tip_racks=[tips_96_10])
	
	# move - 5.0ul from OrderPlate[0] to ReactionPlate[0]
	tips_96_10.pick_up_tip()
	tips_96_10.aspirate(5.0, OrderPlate[0])
	tips_96_10.dispense(5.0, ReactionPlate[0])
	tips_96_10.drop_tip()

	# pause for ['add'] is not currently supported 
	protocol.pause('['add'] is not currently supported ')
	
	# move - 1.0ul from OrderPlate[1] to ReactionPlate[0]
	tips_96_10.pick_up_tip()
	tips_96_10.aspirate(1.0, OrderPlate[1])
	tips_96_10.dispense(1.0, ReactionPlate[0])
	tips_96_10.drop_tip()

	# pause for ['add'] is not currently supported 
	protocol.pause('['add'] is not currently supported ')

	# pause for store product (Cc1ccccc1-c1cccc(NC(=O)Cc2cccc(O)c2)c1)
	protocol.pause('store product (Cc1ccccc1-c1cccc(NC(=O)Cc2cccc(O)c2)c1)')
