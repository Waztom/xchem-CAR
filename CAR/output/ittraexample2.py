from opentrons import protocol_api
# ittraexample2 for "Example Author" produced by XChem Car (https://car.xchem.diamond.ac.uk)
# metadata
 metadata = {
	'protocolName': 'ittraexample2',
	'author': 'Example Author',
	'description': '('N/A',)',
	'apiLevel': '2.9'}

def run(protocol: protocol_api.ProtocolContext):

	# labware
	OrderPlate_1 = protocol.load_labwear('plate_12', '1')
	ReactionPlate_2 = protocol.load_labwear('plate_12', '2')
	tips_96_10_3 = protocol.load_labwear('opentrons_96_filtertiprack_10ul', '3')

	# pipettes
	left_pipette = protocol.load_instrument('p10_single', 'left', tip_racks=[tips_96_10])
	
	# move - 1.0ul from OrderPlate[0] to ReactionPlate[0]
	1.pick_up_tip()
	1.aspirate(1.0, OrderPlate[0])
	1.dispense(1.0, ReactionPlate[0])
	1.drop_tip()

	# pause for ['add'] is not currently supported 
	protocol.pause('['add'] is not currently supported ')
	
	# move - 1.0ul from OrderPlate[1] to ReactionPlate[0]
	1.pick_up_tip()
	1.aspirate(1.0, OrderPlate[1])
	1.dispense(1.0, ReactionPlate[0])
	1.drop_tip()

	# pause for ['add'] is not currently supported 
	protocol.pause('['add'] is not currently supported ')

	# pause for stir at a normal speed at 100.0 Celsius for 3600.0 seconds
	protocol.pause('stir at a normal speed at 100.0 Celsius for 3600.0 seconds')

	# pause for ['set-temperature'] is not currently supported 
	protocol.pause('['set-temperature'] is not currently supported ')

	# pause for store product (Cc1ccccc1-c1cccc(NC(=O)Cc2cccc(O)c2)c1)
	protocol.pause('store product (Cc1ccccc1-c1cccc(NC(=O)Cc2cccc(O)c2)c1)')
