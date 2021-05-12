from opentrons import protocol_api
# ittraexample1 for "Example Author" produced by XChem Car (https://car.xchem.diamond.ac.uk)
# metadata
 metadata = {
	'protocolName': 'ittraexample1',
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
