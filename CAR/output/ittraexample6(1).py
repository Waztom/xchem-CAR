from opentrons import protocol_api
# ittraexample6 for "Example Author" produced by XChem Car (https://car.xchem.diamond.ac.uk)
# metadata
 metadata = {
	'protocolName': 'ittraexample6',
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
	
	# move - 1.0ul from OrderPlate[0] to ReactionPlate[0]
	1.pick_up_tip()
	1.aspirate(1.0, OrderPlate[0])
	1.dispense(1.0, ReactionPlate[0])
	1.drop_tip()

	# pause for ['add'] is not currently supported 
	protocol.pause('['add'] is not currently supported ')
	
	# move - 10.0ul from OrderPlate[1] to ReactionPlate[0]
	1.pick_up_tip()
	1.aspirate(10.0, OrderPlate[1])
	1.dispense(10.0, ReactionPlate[0])
	1.drop_tip()

	# pause for ['add'] is not currently supported 
	protocol.pause('['add'] is not currently supported ')
	
	# move - 1.0ul from OrderPlate[2] to ReactionPlate[0]
	1.pick_up_tip()
	1.aspirate(1.0, OrderPlate[2])
	1.dispense(1.0, ReactionPlate[0])
	1.drop_tip()

	# pause for ['add'] is not currently supported 
	protocol.pause('['add'] is not currently supported ')
	
	# move - 1.0ul from OrderPlate[3] to ReactionPlate[0]
	1.pick_up_tip()
	1.aspirate(1.0, OrderPlate[3])
	1.dispense(1.0, ReactionPlate[0])
	1.drop_tip()

	# pause for ['add'] is not currently supported 
	protocol.pause('['add'] is not currently supported ')
	
	# move - 0.0ul from OrderPlate[] to ReactionPlate[0]
	1.pick_up_tip()
	1.aspirate(0.0, OrderPlate[])
	1.dispense(0.0, ReactionPlate[0])
	1.drop_tip()

	# pause for ['add'] is not currently supported 
	protocol.pause('['add'] is not currently supported ')

	# pause for stir at a normal speed at 25.0 Celsius for 86400.0 seconds
	protocol.pause('stir at a normal speed at 25.0 Celsius for 86400.0 seconds')
	
	# move - 10.0ul from OrderPlate[4] to ReactionPlate[0]
	1.pick_up_tip()
	1.aspirate(10.0, OrderPlate[4])
	1.dispense(10.0, ReactionPlate[0])
	1.drop_tip()

	# pause for ['add'] is not currently supported 
	protocol.pause('['add'] is not currently supported ')

	# pause for extract (nan)
	protocol.pause('extract (nan)')

	# pause for extract (nan)
	protocol.pause('extract (nan)')
	
	# move - nanul from OrderPlate[] to ReactionPlate[0]
	1.pick_up_tip()
	1.aspirate(nan, OrderPlate[])
	1.dispense(nan, ReactionPlate[0])
	1.drop_tip()

	# pause for ['collect-layer'] is not currently supported 
	protocol.pause('['collect-layer'] is not currently supported ')
	
	# move - 10.0ul from OrderPlate[] to ReactionPlate[0]
	1.pick_up_tip()
	1.aspirate(10.0, OrderPlate[])
	1.dispense(10.0, ReactionPlate[0])
	1.drop_tip()

	# pause for [] is not currently supported 
	protocol.pause('[] is not currently supported ')

	# pause for concentrate (nan)
	protocol.pause('concentrate (nan)')
