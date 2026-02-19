"""Tests for backend.db_utils module.

Covers:
  - Action QuerySets (AddAction, ExtractAction, OTBatchProtocol)
  - ActionSession utilities (sequence numbers, types, grouping, query sets)
  - Query helpers (previous object entries)
  - Reaction checks (failure checks, method steps, reactions to do)
  - Basic model queries (targets, methods, reactions, batch operations)
  - Batch data retrieval (SMILES, reaction IDs, product SMILES)
  - Update operations (target mols, reaction success, OT-friendly, recipe, etc.)
  - Plate queries (query sets, plate maps)
  - Product queries (query sets, single product, product SMILES)
  - Reaction queries (single, temperature, class, recipe, query sets, proceeding)
  - Reaction yields and previous reaction products
  - Batch MW calculations (target MWs, reaction product MWs)

All external dependencies (Django ORM, chem_utils, recipe_utils, conversions)
are mocked.
"""

from unittest import TestCase
from unittest.mock import patch, MagicMock, call, mock_open

from backend.db_utils import (
    get_add_action_query_set,
    get_extract_action_query_set,
    get_ot_batch_protocol_query_set,
    get_action_session_sequence_numbers,
    get_action_session_types,
    get_grouped_action_session_sequences,
    get_grouped_action_session_types,
    get_action_session_query_set,
    get_previous_obj_entries,
    check_previous_reaction_failures,
    check_no_method_steps,
    get_reactions_to_do,
    get_targets,
    get_methods,
    get_reactions,
    get_batch_tag,
    get_batch_reactions,
    get_max_reaction_number,
    group_reactions,
    get_reactants_to_buy,
    get_batch_target_smiles,
    get_batch_reaction_ids,
    get_batch_reaction_product_smiles,
    update_target_mols,
    update_reaction_success_to_fail,
    update_batch_method_ot_friendly,
    delete_batch_action_sessions,
    update_recipe_type,
    get_plate_query_set,
    get_plate_map,
    get_product_query_set,
    get_product,
    get_product_smiles,
    get_reaction,
    get_reaction_temperature,
    get_reaction_class,
    get_reaction_recipe,
    get_reaction_query_set,
    check_proceeding_reactions,
    get_reaction_yields,
    check_previous_reaction_products,
    get_previous_reaction_query_sets,
    get_batch_target_mws,
    get_batch_reaction_product_mws,
)

# ---------------------------------------------------------------------------
# Module path constant used for patching
# ---------------------------------------------------------------------------
_DB = "backend.db_utils"


# ===========================================================================
# Action QuerySets
# ===========================================================================


class TestGetAddActionQuerySet(TestCase):
    """Tests for get_add_action_query_set."""

    @patch(f"{_DB}.AddAction")
    def test_with_actionsession_ids(self, MockAddAction):
        mock_qs = MagicMock()
        MockAddAction.objects.filter.return_value.order_by.return_value = mock_qs

        result = get_add_action_query_set(
            reaction_ids=[1, 2], actionsession_ids=[10, 20]
        )

        self.assertIs(result, mock_qs)
        MockAddAction.objects.filter.assert_called_once()

    @patch(f"{_DB}.AddAction")
    def test_with_actionsessiontype(self, MockAddAction):
        mock_qs = MagicMock()
        MockAddAction.objects.filter.return_value.order_by.return_value = mock_qs

        result = get_add_action_query_set(
            reaction_ids=[1], actionsessiontype="reaction"
        )

        self.assertIs(result, mock_qs)
        MockAddAction.objects.filter.assert_called_once()

    @patch(f"{_DB}.AddAction")
    def test_no_optional_args_returns_none(self, MockAddAction):
        result = get_add_action_query_set(reaction_ids=[1, 2])
        self.assertIsNone(result)
        MockAddAction.objects.filter.assert_not_called()


class TestGetExtractActionQuerySet(TestCase):
    """Tests for get_extract_action_query_set."""

    @patch(f"{_DB}.ExtractAction")
    def test_with_actionsession_ids(self, MockExtractAction):
        mock_qs = MagicMock()
        MockExtractAction.objects.filter.return_value.order_by.return_value = mock_qs

        result = get_extract_action_query_set(
            reaction_ids=[3, 4], actionsession_ids=[30]
        )

        self.assertIs(result, mock_qs)

    @patch(f"{_DB}.ExtractAction")
    def test_with_actionsessiontype(self, MockExtractAction):
        mock_qs = MagicMock()
        MockExtractAction.objects.filter.return_value.order_by.return_value = mock_qs

        result = get_extract_action_query_set(
            reaction_ids=[5], actionsessiontype="workup"
        )

        self.assertIs(result, mock_qs)

    @patch(f"{_DB}.ExtractAction")
    def test_no_optional_args_returns_none(self, MockExtractAction):
        result = get_extract_action_query_set(reaction_ids=[1])
        self.assertIsNone(result)


class TestGetOTBatchProtocolQuerySet(TestCase):
    """Tests for get_ot_batch_protocol_query_set."""

    @patch(f"{_DB}.OTBatchProtocol")
    def test_returns_filtered_queryset(self, MockOTBP):
        mock_qs = MagicMock()
        MockOTBP.objects.filter.return_value = mock_qs

        result = get_ot_batch_protocol_query_set(batch_id=42)

        self.assertIs(result, mock_qs)
        MockOTBP.objects.filter.assert_called_once_with(batch_id=42)


# ===========================================================================
# ActionSession utilities
# ===========================================================================


class TestGetActionSessionSequenceNumbers(TestCase):
    """Tests for get_action_session_sequence_numbers."""

    def test_returns_range_list(self):
        qs = MagicMock()
        qs.aggregate.return_value = {"sessionnumber__max": 4}

        result = get_action_session_sequence_numbers(qs)

        self.assertEqual(result, [1, 2, 3, 4])
        qs.aggregate.assert_called_once()

    def test_single_session(self):
        qs = MagicMock()
        qs.aggregate.return_value = {"sessionnumber__max": 1}

        result = get_action_session_sequence_numbers(qs)
        self.assertEqual(result, [1])


class TestGetActionSessionTypes(TestCase):
    """Tests for get_action_session_types."""

    def test_returns_set_of_types(self):
        qs = MagicMock()
        qs.values_list.return_value = ["reaction", "workup", "reaction"]

        result = get_action_session_types(qs)

        self.assertEqual(result, {"reaction", "workup"})
        qs.values_list.assert_called_once_with("type", flat=True)

    def test_single_type(self):
        qs = MagicMock()
        qs.values_list.return_value = ["stir"]

        result = get_action_session_types(qs)
        self.assertEqual(result, {"stir"})


class TestGetGroupedActionSessionSequences(TestCase):
    """Tests for get_grouped_action_session_sequences."""

    def test_groups_by_sequence_number(self):
        qs = MagicMock()
        group1 = MagicMock(name="group1")
        group2 = MagicMock(name="group2")
        qs.filter.return_value.order_by.side_effect = [group1, group2]

        result = get_grouped_action_session_sequences([1, 2], qs)

        self.assertEqual(len(result), 2)
        self.assertIs(result[0], group1)
        self.assertIs(result[1], group2)

    def test_empty_session_numbers(self):
        qs = MagicMock()
        result = get_grouped_action_session_sequences([], qs)
        self.assertEqual(result, [])


class TestGetGroupedActionSessionTypes(TestCase):
    """Tests for get_grouped_action_session_types."""

    def test_groups_by_type(self):
        qs = MagicMock()
        group1 = MagicMock(__bool__=lambda s: True)
        group2 = MagicMock(__bool__=lambda s: True)
        qs.filter.return_value.order_by.side_effect = [group1, group2]

        result = get_grouped_action_session_types({"reaction", "workup"}, qs)

        self.assertEqual(len(result), 2)

    def test_skips_empty_groups(self):
        qs = MagicMock()
        empty = MagicMock(__bool__=lambda s: False)
        filled = MagicMock(__bool__=lambda s: True)
        qs.filter.return_value.order_by.side_effect = [empty, filled]

        result = get_grouped_action_session_types({"reaction", "workup"}, qs)

        self.assertEqual(len(result), 1)
        self.assertIs(result[0], filled)


class TestGetActionSessionQuerySet(TestCase):
    """Tests for get_action_session_query_set."""

    @patch(f"{_DB}.ActionSession")
    def test_with_driver_found(self, MockAS):
        mock_qs = MagicMock(__bool__=lambda s: True)
        MockAS.objects.filter.return_value = mock_qs

        result = get_action_session_query_set(
            reaction_ids=[1, 2], driver="robot"
        )

        self.assertIs(result, mock_qs)

    @patch(f"{_DB}.ActionSession")
    def test_without_driver(self, MockAS):
        mock_qs = MagicMock(__bool__=lambda s: True)
        MockAS.objects.filter.return_value = mock_qs

        result = get_action_session_query_set(reaction_ids=[1, 2])

        self.assertIs(result, mock_qs)

    @patch(f"{_DB}.ActionSession")
    def test_with_driver_not_found_returns_none(self, MockAS):
        """When driver is given but queryset is empty, returns None."""
        MockAS.objects.filter.return_value = MagicMock(__bool__=lambda s: False)

        result = get_action_session_query_set(
            reaction_ids=[1], driver="robot"
        )

        self.assertIsNone(result)

    @patch(f"{_DB}.ActionSession")
    def test_no_driver_no_results_returns_none(self, MockAS):
        MockAS.objects.filter.return_value = MagicMock(__bool__=lambda s: False)

        result = get_action_session_query_set(reaction_ids=[1])
        self.assertIsNone(result)


# ===========================================================================
# Query helpers
# ===========================================================================


class TestGetPreviousObjEntries(TestCase):
    """Tests for get_previous_obj_entries."""

    def test_filters_and_orders(self):
        qs = MagicMock()
        obj = MagicMock(pk=10)
        expected = MagicMock()
        qs.filter.return_value.order_by.return_value = expected

        result = get_previous_obj_entries(qs, obj)

        self.assertIs(result, expected)
        qs.filter.assert_called_once_with(pk__lt=10)
        qs.filter.return_value.order_by.assert_called_once_with("-pk")


# ===========================================================================
# Reaction checks and queries
# ===========================================================================


class TestCheckPreviousReactionFailures(TestCase):
    """Tests for check_previous_reaction_failures."""

    @patch(f"{_DB}.get_previous_obj_entries")
    @patch(f"{_DB}.get_reactions")
    def test_returns_true_when_failures_exist(self, mock_get_rxn, mock_prev):
        rxn = MagicMock()
        rxn.method_id.id = 5

        mock_get_rxn.return_value = MagicMock()
        prev_qs = MagicMock()
        mock_prev.return_value = prev_qs
        prev_qs.filter.return_value.exists.return_value = True

        self.assertTrue(check_previous_reaction_failures(rxn))

    @patch(f"{_DB}.get_previous_obj_entries")
    @patch(f"{_DB}.get_reactions")
    def test_returns_false_when_no_failures(self, mock_get_rxn, mock_prev):
        rxn = MagicMock()
        rxn.method_id.id = 5

        mock_get_rxn.return_value = MagicMock()
        prev_qs = MagicMock()
        mock_prev.return_value = prev_qs
        prev_qs.filter.return_value.exists.return_value = False

        self.assertFalse(check_previous_reaction_failures(rxn))


class TestCheckNoMethodSteps(TestCase):
    """Tests for check_no_method_steps."""

    def test_returns_true_when_more_than_one_step(self):
        rxn = MagicMock()
        rxn.method_id.nosteps = 3
        self.assertTrue(check_no_method_steps(rxn))

    def test_returns_false_when_single_step(self):
        rxn = MagicMock()
        rxn.method_id.nosteps = 1
        self.assertFalse(check_no_method_steps(rxn))


class TestGetReactionsToDo(TestCase):
    """Tests for get_reactions_to_do."""

    @patch(f"{_DB}.check_previous_reaction_failures")
    @patch(f"{_DB}.check_no_method_steps")
    def test_excludes_failed_reactions(self, mock_steps, mock_fail):
        rxn1 = MagicMock(id=1)
        rxn2 = MagicMock(id=2)
        rxn3 = MagicMock(id=3)
        group_qs = MagicMock()
        group_qs.__iter__ = lambda s: iter([rxn1, rxn2, rxn3])

        # rxn1: multi-step + no failures -> include
        # rxn2: multi-step + previous failure -> exclude
        # rxn3: single-step -> exclude
        mock_steps.side_effect = [True, True, False]
        mock_fail.side_effect = [False, True]

        filtered_qs = MagicMock()
        group_qs.filter.return_value = filtered_qs

        result = get_reactions_to_do(group_qs)

        group_qs.filter.assert_called_once_with(id__in=[1])
        self.assertIs(result, filtered_qs)

    @patch(f"{_DB}.check_previous_reaction_failures")
    @patch(f"{_DB}.check_no_method_steps")
    def test_all_single_step_returns_empty(self, mock_steps, mock_fail):
        rxn = MagicMock(id=1)
        group_qs = MagicMock()
        group_qs.__iter__ = lambda s: iter([rxn])
        mock_steps.return_value = False

        filtered_qs = MagicMock()
        group_qs.filter.return_value = filtered_qs

        result = get_reactions_to_do(group_qs)
        group_qs.filter.assert_called_once_with(id__in=[])


# ===========================================================================
# Basic model queries
# ===========================================================================


class TestGetTargets(TestCase):
    """Tests for get_targets."""

    @patch(f"{_DB}.Target")
    def test_returns_ordered_queryset(self, MockTarget):
        mock_qs = MagicMock()
        MockTarget.objects.filter.return_value.order_by.return_value = mock_qs

        result = get_targets(batch_ids=[1, 2])

        self.assertIs(result, mock_qs)
        MockTarget.objects.filter.assert_called_once_with(batch_id__in=[1, 2])


class TestGetMethods(TestCase):
    """Tests for get_methods."""

    @patch(f"{_DB}.Method")
    def test_filters_otchem_and_orders(self, MockMethod):
        mock_qs = MagicMock()
        (
            MockMethod.objects.filter.return_value
            .filter.return_value
            .order_by.return_value
        ) = mock_qs

        result = get_methods(target_ids=[10])

        self.assertIs(result, mock_qs)


class TestGetReactions(TestCase):
    """Tests for get_reactions."""

    @patch(f"{_DB}.Reaction")
    def test_returns_ordered_queryset(self, MockReaction):
        mock_qs = MagicMock()
        MockReaction.objects.filter.return_value.order_by.return_value = mock_qs

        result = get_reactions(method_ids=[5])

        self.assertIs(result, mock_qs)
        MockReaction.objects.filter.assert_called_once_with(method_id__in=[5])


class TestGetBatchTag(TestCase):
    """Tests for get_batch_tag."""

    @patch(f"{_DB}.Batch")
    def test_returns_batch_tag(self, MockBatch):
        batch_obj = MagicMock()
        batch_obj.batchtag = "BATCH-001"
        MockBatch.objects.get.return_value = batch_obj

        result = get_batch_tag(batchid=1)

        self.assertEqual(result, "BATCH-001")
        MockBatch.objects.get.assert_called_once_with(id=1)


class TestGetBatchReactions(TestCase):
    """Tests for get_batch_reactions."""

    @patch(f"{_DB}.get_reactions")
    @patch(f"{_DB}.get_methods")
    @patch(f"{_DB}.get_targets")
    def test_happy_path(self, mock_targets, mock_methods, mock_rxns):
        mock_targets.return_value = MagicMock(__bool__=lambda s: True)
        mock_methods.return_value = MagicMock(__bool__=lambda s: True)
        expected = MagicMock(__bool__=lambda s: True)
        mock_rxns.return_value = expected

        result = get_batch_reactions(batchid=1)

        self.assertIs(result, expected)

    @patch(f"{_DB}.get_targets")
    def test_no_targets_returns_none(self, mock_targets):
        mock_targets.return_value = MagicMock(__bool__=lambda s: False)

        result = get_batch_reactions(batchid=1)

        self.assertIsNone(result)

    @patch(f"{_DB}.get_methods")
    @patch(f"{_DB}.get_targets")
    def test_no_methods_returns_none(self, mock_targets, mock_methods):
        mock_targets.return_value = MagicMock(__bool__=lambda s: True)
        mock_methods.return_value = MagicMock(__bool__=lambda s: False)

        result = get_batch_reactions(batchid=1)

        self.assertIsNone(result)

    @patch(f"{_DB}.get_reactions")
    @patch(f"{_DB}.get_methods")
    @patch(f"{_DB}.get_targets")
    def test_no_reactions_returns_none(self, mock_t, mock_m, mock_r):
        mock_t.return_value = MagicMock(__bool__=lambda s: True)
        mock_m.return_value = MagicMock(__bool__=lambda s: True)
        mock_r.return_value = MagicMock(__bool__=lambda s: False)

        result = get_batch_reactions(batchid=1)

        self.assertIsNone(result)


class TestGetMaxReactionNumber(TestCase):
    """Tests for get_max_reaction_number."""

    def test_returns_max(self):
        qs = MagicMock()
        qs.aggregate.return_value = {"number__max": 3}

        result = get_max_reaction_number(qs)

        self.assertEqual(result, 3)


class TestGroupReactions(TestCase):
    """Tests for group_reactions."""

    def test_groups_by_number(self):
        qs = MagicMock()
        group1 = MagicMock(__bool__=lambda s: True)
        group2 = MagicMock(__bool__=lambda s: True)
        qs.filter.return_value.distinct.return_value.order_by.side_effect = [
            group1, group2,
        ]

        result = group_reactions(qs, maxreactionnumber=2)

        self.assertEqual(len(result), 2)
        self.assertIs(result[0], group1)
        self.assertIs(result[1], group2)

    def test_skips_empty_groups(self):
        qs = MagicMock()
        filled = MagicMock(__bool__=lambda s: True)
        empty = MagicMock(__bool__=lambda s: False)
        qs.filter.return_value.distinct.return_value.order_by.side_effect = [
            empty, filled,
        ]

        result = group_reactions(qs, maxreactionnumber=2)

        self.assertEqual(len(result), 1)
        self.assertIs(result[0], filled)


class TestGetReactantsToBuy(TestCase):
    """Tests for get_reactants_to_buy."""

    @patch(f"{_DB}.Batch")
    def test_returns_unique_smiles(self, MockBatch):
        reactant1 = MagicMock(smiles="CCO", previousreactionproduct=False)
        reactant2 = MagicMock(smiles="CC=O", previousreactionproduct=False)
        reactant3 = MagicMock(smiles="CCO", previousreactionproduct=True)

        rxn = MagicMock()
        rxn.reactants.all.return_value = [reactant1, reactant2, reactant3]

        method = MagicMock()
        method.reactions.all.return_value = [rxn]

        target = MagicMock()
        target.methods.all.return_value = [method]

        batch = MagicMock()
        batch.targets.all.return_value = [target]
        MockBatch.objects.get.return_value = batch

        result = get_reactants_to_buy(batch_ids=[1])

        self.assertIsInstance(result, list)
        self.assertEqual(set(result), {"CCO", "CC=O"})

    @patch(f"{_DB}.Batch")
    def test_all_previous_reaction_products(self, MockBatch):
        """When all reactants are previous-reaction products, list is empty."""
        reactant = MagicMock(smiles="CCO", previousreactionproduct=True)
        rxn = MagicMock()
        rxn.reactants.all.return_value = [reactant]
        method = MagicMock()
        method.reactions.all.return_value = [rxn]
        target = MagicMock()
        target.methods.all.return_value = [method]
        batch = MagicMock()
        batch.targets.all.return_value = [target]
        MockBatch.objects.get.return_value = batch

        result = get_reactants_to_buy(batch_ids=[1])
        self.assertEqual(result, [])


# ===========================================================================
# Batch data retrieval
# ===========================================================================


class TestGetBatchTargetSmiles(TestCase):
    """Tests for get_batch_target_smiles."""

    @patch(f"{_DB}.Batch")
    def test_returns_smiles_list(self, MockBatch):
        t1 = MagicMock(smiles="CCO")
        t2 = MagicMock(smiles="CC=O")
        batch = MagicMock()
        batch.targets.all.return_value = [t1, t2]
        MockBatch.objects.get.return_value = batch

        result = get_batch_target_smiles(batch_id=1)

        self.assertEqual(result, ["CCO", "CC=O"])


class TestGetBatchReactionIds(TestCase):
    """Tests for get_batch_reaction_ids."""

    @patch(f"{_DB}.Batch")
    def test_returns_ids(self, MockBatch):
        r1 = MagicMock(id=100)
        r2 = MagicMock(id=200)

        method = MagicMock()
        method.reactions.all.return_value.filter.return_value.order_by.return_value = [
            r1, r2,
        ]

        target = MagicMock()
        target.methods.all.return_value.order_by.return_value = [method]

        batch = MagicMock()
        batch.targets.all.return_value.order_by.return_value = [target]
        MockBatch.objects.get.return_value = batch

        result = get_batch_reaction_ids(batch_id=1, reaction_number=1)

        self.assertEqual(result, [100, 200])


class TestGetBatchReactionProductSmiles(TestCase):
    """Tests for get_batch_reaction_product_smiles."""

    @patch(f"{_DB}.Batch")
    def test_returns_product_smiles(self, MockBatch):
        rxn = MagicMock()
        rxn.products.all.return_value.values_list.return_value = ["CCO", "CC=O"]

        method = MagicMock()
        method.reactions.all.return_value.filter.return_value.order_by.return_value = [
            rxn
        ]

        target = MagicMock()
        target.methods.all.return_value.order_by.return_value = [method]

        batch = MagicMock()
        batch.targets.all.return_value.order_by.return_value = [target]
        MockBatch.objects.get.return_value = batch

        result = get_batch_reaction_product_smiles(batch_id=1, reaction_number=1)

        self.assertEqual(result, ["CCO", "CC=O"])


# ===========================================================================
# Update operations
# ===========================================================================


class TestUpdateTargetMols(TestCase):
    """Tests for update_target_mols."""

    @patch(f"{_DB}.Batch")
    @patch(f"{_DB}.calculate_mass_from_mols", create=True)
    @patch(f"{_DB}.calculate_mols_from_conc", create=True)
    def test_updates_targets(self, mock_mols_calc, mock_mass_calc, MockBatch):
        """Patch at the deferred-import location inside the function."""
        target = MagicMock(smiles="CCO")
        MockBatch.objects.get.return_value.targets.all.return_value = [target]

        with patch(
            "backend.conversions.calculate_mols_from_conc", return_value=0.001
        ) as mock_mc, patch(
            "backend.conversions.calculate_mass_from_mols", return_value=0.046
        ) as mock_mm:
            update_target_mols(batch_id=1, concentration=10.0, volume=100.0)

        target.save.assert_called_once()


class TestUpdateReactionSuccessToFail(TestCase):
    """Tests for update_reaction_success_to_fail."""

    @patch(f"{_DB}.Reaction")
    def test_updates_when_exists(self, MockReaction):
        MockReaction.objects.filter.return_value.exists.return_value = True

        update_reaction_success_to_fail(reaction_ids=[1, 2])

        # Called twice: once for exists(), once for update()
        self.assertEqual(MockReaction.objects.filter.call_count, 2)

    @patch(f"{_DB}.Reaction")
    def test_noop_when_not_exists(self, MockReaction):
        MockReaction.objects.filter.return_value.exists.return_value = False

        update_reaction_success_to_fail(reaction_ids=[99])

        MockReaction.objects.filter.return_value.update.assert_not_called()


class TestUpdateBatchMethodOtFriendly(TestCase):
    """Tests for update_batch_method_ot_friendly."""

    @patch(f"{_DB}.Batch")
    def test_sets_otchem_true(self, MockBatch):
        method = MagicMock()
        target = MagicMock()
        target.methods.all.return_value.order_by.return_value = [method]
        batch = MagicMock()
        batch.targets.all.return_value.order_by.return_value = [target]
        MockBatch.objects.get.return_value = batch

        update_batch_method_ot_friendly(batch_id=1)

        self.assertTrue(method.otchem)
        method.save.assert_called_once()


class TestDeleteBatchActionSessions(TestCase):
    """Tests for delete_batch_action_sessions."""

    @patch(f"{_DB}.get_action_session_query_set")
    @patch(f"{_DB}.get_batch_reactions")
    def test_deletes_sessions(self, mock_batch_rxns, mock_as_qs):
        rxn_qs = MagicMock()
        mock_batch_rxns.return_value = rxn_qs
        as_qs = MagicMock()
        mock_as_qs.return_value = as_qs

        delete_batch_action_sessions(batch_id=1)

        mock_batch_rxns.assert_called_once_with(batchid=1)
        mock_as_qs.assert_called_once_with(reaction_ids=rxn_qs)
        as_qs.delete.assert_called_once()


class TestUpdateRecipeType(TestCase):
    """Tests for update_recipe_type."""

    @patch(f"{_DB}.get_batch_reactions")
    def test_updates_recipes(self, mock_batch_rxns):
        rxn1 = MagicMock(recipe="old_recipe")
        rxn2 = MagicMock(recipe="old_recipe")
        qs = MagicMock()
        qs.__iter__ = lambda s: iter([rxn1, rxn2])
        mock_batch_rxns.return_value.filter.return_value = qs

        update_recipe_type(
            batch_id=1,
            reaction_class="amidation",
            current_recipe="old_recipe",
            recipe_to_use="new_recipe",
        )

        self.assertEqual(rxn1.recipe, "new_recipe")
        self.assertEqual(rxn2.recipe, "new_recipe")
        rxn1.save.assert_called_once()
        rxn2.save.assert_called_once()


# ===========================================================================
# Plate queries
# ===========================================================================


class TestGetPlateQuerySet(TestCase):
    """Tests for get_plate_query_set."""

    @patch(f"{_DB}.Plate")
    def test_by_plate_id(self, MockPlate):
        mock_qs = MagicMock()
        MockPlate.objects.filter.return_value = mock_qs

        result = get_plate_query_set(plate_id=5)

        self.assertIs(result, mock_qs)
        MockPlate.objects.filter.assert_called_with(id=5)

    @patch(f"{_DB}.Plate")
    def test_by_otsession_id(self, MockPlate):
        mock_qs = MagicMock()
        MockPlate.objects.filter.return_value = mock_qs

        result = get_plate_query_set(otsession_id=10)

        self.assertIs(result, mock_qs)
        MockPlate.objects.filter.assert_called_with(otsession_id=10)


class TestGetPlateMap(TestCase):
    """Tests for get_plate_map."""

    @patch("builtins.open", new_callable=mock_open)
    @patch(f"{_DB}.csv")
    @patch(f"{_DB}.Plate")
    def test_writes_csv(self, MockPlate, mock_csv, mock_file):
        well = MagicMock()
        well.index = "A1"
        well.method_id.target_id.id = 1
        well.method_id.target_id.name = "target1"
        well.smiles = "CCO"
        well.reaction_id.reactants.values_list.return_value = ["CCO", "CC=O"]

        plate = MagicMock()
        plate.well_set.all.return_value.order_by.return_value = [well]
        MockPlate.objects.get.return_value = plate

        with patch(
            "backend.chem_utils.get_mws", return_value=[46.07, 44.05]
        ):
            get_plate_map(plate_ids=[1], out_dir="/tmp/")

        mock_csv.writer.assert_called_once()
        writer = mock_csv.writer.return_value
        writer.writerow.assert_called()

    @patch(f"{_DB}.Plate")
    def test_exception_logs_error(self, MockPlate):
        MockPlate.objects.get.side_effect = Exception("not found")

        # Should not raise — exception is caught internally
        get_plate_map(plate_ids=[99], out_dir="/tmp/")


# ===========================================================================
# Product queries
# ===========================================================================


class TestGetProductQuerySet(TestCase):
    """Tests for get_product_query_set."""

    @patch(f"{_DB}.Product")
    def test_returns_filtered_queryset(self, MockProduct):
        mock_qs = MagicMock()
        MockProduct.objects.filter.return_value = mock_qs

        result = get_product_query_set(reaction_ids=[1, 2])

        self.assertIs(result, mock_qs)
        MockProduct.objects.filter.assert_called_once_with(
            reaction_id__in=[1, 2]
        )


class TestGetProduct(TestCase):
    """Tests for get_product."""

    @patch(f"{_DB}.Product")
    def test_returns_product_object(self, MockProduct):
        product = MagicMock()
        MockProduct.objects.get.return_value = product

        result = get_product(reaction_id=5)

        self.assertIs(result, product)
        MockProduct.objects.get.assert_called_once_with(reaction_id=5)


class TestGetProductSmiles(TestCase):
    """Tests for get_product_smiles."""

    @patch(f"{_DB}.Product")
    def test_returns_list_when_found(self, MockProduct):
        vl = MagicMock(__bool__=lambda s: True, __iter__=lambda s: iter(["CCO"]))
        MockProduct.objects.filter.return_value.values_list.return_value = vl

        result = get_product_smiles(reaction_ids=[1])

        self.assertEqual(result, ["CCO"])

    @patch(f"{_DB}.Product")
    def test_returns_none_when_empty(self, MockProduct):
        vl = MagicMock(__bool__=lambda s: False)
        MockProduct.objects.filter.return_value.values_list.return_value = vl

        result = get_product_smiles(reaction_ids=[999])

        self.assertIsNone(result)


# ===========================================================================
# Reaction queries
# ===========================================================================


class TestGetReaction(TestCase):
    """Tests for get_reaction."""

    @patch(f"{_DB}.Reaction")
    def test_returns_reaction(self, MockReaction):
        rxn = MagicMock()
        MockReaction.objects.get.return_value = rxn

        result = get_reaction(reaction_id=7)

        self.assertIs(result, rxn)
        MockReaction.objects.get.assert_called_once_with(id=7)


class TestGetReactionTemperature(TestCase):
    """Tests for get_reaction_temperature."""

    @patch(f"{_DB}.Reaction")
    def test_returns_temperature(self, MockReaction):
        rxn = MagicMock(temperature=80.0)
        MockReaction.objects.get.return_value = rxn

        result = get_reaction_temperature(reaction_id=7)

        self.assertEqual(result, 80.0)


class TestGetReactionClass(TestCase):
    """Tests for get_reaction_class."""

    @patch(f"{_DB}.Reaction")
    def test_returns_class(self, MockReaction):
        rxn = MagicMock(reactionclass="amidation")
        MockReaction.objects.get.return_value = rxn

        result = get_reaction_class(reaction_id=7)

        self.assertEqual(result, "amidation")


class TestGetReactionRecipe(TestCase):
    """Tests for get_reaction_recipe."""

    @patch(f"{_DB}.Reaction")
    def test_returns_recipe(self, MockReaction):
        rxn = MagicMock(recipe="recipe_v2")
        MockReaction.objects.get.return_value = rxn

        result = get_reaction_recipe(reaction_id=7)

        self.assertEqual(result, "recipe_v2")


class TestGetReactionQuerySet(TestCase):
    """Tests for get_reaction_query_set."""

    @patch(f"{_DB}.Reaction")
    def test_by_reaction_ids(self, MockReaction):
        mock_qs = MagicMock()
        MockReaction.objects.filter.return_value.order_by.return_value = mock_qs

        result = get_reaction_query_set(reaction_ids=[1, 2])

        self.assertIs(result, mock_qs)

    @patch(f"{_DB}.Reaction")
    def test_by_method_id(self, MockReaction):
        mock_qs = MagicMock()
        MockReaction.objects.filter.return_value.order_by.return_value = mock_qs

        result = get_reaction_query_set(method_id=10)

        self.assertIs(result, mock_qs)
        MockReaction.objects.filter.assert_called_once_with(method_id=10)


class TestCheckProceedingReactions(TestCase):
    """Tests for check_proceeding_reactions."""

    @patch(f"{_DB}.Method")
    @patch(f"{_DB}.get_reaction")
    def test_returns_proceeding_queryset(self, mock_get_rxn, MockMethod):
        rxn = MagicMock()
        rxn.method_id.id = 5
        mock_get_rxn.return_value = rxn

        proceeding_qs = MagicMock()
        MockMethod.objects.get.return_value.reactions.filter.return_value = (
            proceeding_qs
        )

        result = check_proceeding_reactions(reaction_id=10)

        self.assertIs(result, proceeding_qs)
        mock_get_rxn.assert_called_once_with(reaction_id=10)
        MockMethod.objects.get.assert_called_once_with(id=5)
        MockMethod.objects.get.return_value.reactions.filter.assert_called_once_with(
            id__gt=10
        )


class TestGetReactionYields(TestCase):
    """Tests for get_reaction_yields."""

    def test_returns_yields(self):
        with patch(
            "backend.recipe_utils.get_recipe_yield", side_effect=[0.8, 0.6]
        ):
            result = get_reaction_yields(
                reactionclasslist=["amidation", "SNAr"],
                recipelist=["recipe1", "recipe2"],
            )

        self.assertEqual(result, [0.8, 0.6])


class TestCheckPreviousReactionProducts(TestCase):
    """Tests for check_previous_reaction_products."""

    @patch(f"{_DB}.get_product")
    @patch(f"{_DB}.get_previous_obj_entries")
    @patch(f"{_DB}.get_reaction_query_set")
    @patch(f"{_DB}.get_reaction")
    def test_returns_true_when_match(
        self, mock_get_rxn, mock_rxn_qs, mock_prev, mock_product
    ):
        rxn = MagicMock()
        rxn.method_id.id = 5
        mock_get_rxn.return_value = rxn

        mock_rxn_qs.return_value = MagicMock()

        prev_rxn = MagicMock()
        prev_qs = MagicMock(
            __bool__=lambda s: True,
            __iter__=lambda s: iter([prev_rxn]),
        )
        mock_prev.return_value = prev_qs

        product = MagicMock(smiles="CCO")
        mock_product.return_value = product

        result = check_previous_reaction_products(reaction_id=10, smiles="CCO")

        self.assertTrue(result)

    @patch(f"{_DB}.get_product")
    @patch(f"{_DB}.get_previous_obj_entries")
    @patch(f"{_DB}.get_reaction_query_set")
    @patch(f"{_DB}.get_reaction")
    def test_returns_false_when_no_match(
        self, mock_get_rxn, mock_rxn_qs, mock_prev, mock_product
    ):
        rxn = MagicMock()
        rxn.method_id.id = 5
        mock_get_rxn.return_value = rxn
        mock_rxn_qs.return_value = MagicMock()

        prev_rxn = MagicMock()
        prev_qs = MagicMock(
            __bool__=lambda s: True,
            __iter__=lambda s: iter([prev_rxn]),
        )
        mock_prev.return_value = prev_qs

        product = MagicMock(smiles="DIFFERENT")
        mock_product.return_value = product

        result = check_previous_reaction_products(reaction_id=10, smiles="CCO")

        self.assertFalse(result)

    @patch(f"{_DB}.get_previous_obj_entries")
    @patch(f"{_DB}.get_reaction_query_set")
    @patch(f"{_DB}.get_reaction")
    def test_returns_false_when_no_previous(
        self, mock_get_rxn, mock_rxn_qs, mock_prev
    ):
        rxn = MagicMock()
        rxn.method_id.id = 5
        mock_get_rxn.return_value = rxn
        mock_rxn_qs.return_value = MagicMock()
        mock_prev.return_value = MagicMock(__bool__=lambda s: False)

        result = check_previous_reaction_products(reaction_id=10, smiles="CCO")

        self.assertFalse(result)


class TestGetPreviousReactionQuerySets(TestCase):
    """Tests for get_previous_reaction_query_sets."""

    @patch(f"{_DB}.Method")
    @patch(f"{_DB}.get_reaction")
    def test_returns_previous_matching(self, mock_get_rxn, MockMethod):
        rxn = MagicMock()
        rxn.method_id.id = 5
        mock_get_rxn.return_value = rxn

        prev_qs = MagicMock()
        MockMethod.objects.get.return_value.reactions.filter.return_value = prev_qs

        result = get_previous_reaction_query_sets(reaction_id=10, smiles="CCO")

        self.assertIs(result, prev_qs)
        MockMethod.objects.get.return_value.reactions.filter.assert_called_once_with(
            id__lt=10, products__smiles="CCO"
        )


# ===========================================================================
# Batch MW calculations
# ===========================================================================


class TestGetBatchTargetMws(TestCase):
    """Tests for get_batch_target_mws."""

    @patch(f"{_DB}.Batch")
    def test_returns_mws(self, MockBatch):
        t1 = MagicMock(smiles="CCO")
        t2 = MagicMock(smiles="CC=O")
        batch = MagicMock()
        batch.targets.all.return_value = [t1, t2]
        MockBatch.objects.get.return_value = batch

        with patch(
            "backend.chem_utils.get_mws", return_value=[46.07, 44.05]
        ):
            result = get_batch_target_mws(batch_id=1)

        self.assertEqual(result, [46.07, 44.05])


class TestGetBatchReactionProductMws(TestCase):
    """Tests for get_batch_reaction_product_mws."""

    @patch(f"{_DB}.Batch")
    def test_returns_product_mws(self, MockBatch):
        rxn = MagicMock()
        rxn.products.all.return_value.order_by.return_value.values_list.return_value = [
            "CCO"
        ]

        method = MagicMock()
        method.reactions.all.return_value.filter.return_value.order_by.return_value = [
            rxn
        ]

        target = MagicMock()
        target.methods.all.return_value.order_by.return_value = [method]

        batch = MagicMock()
        batch.targets.all.return_value.order_by.return_value = [target]
        MockBatch.objects.get.return_value = batch

        with patch(
            "backend.chem_utils.get_mws", return_value=[46.07]
        ):
            result = get_batch_reaction_product_mws(batch_id=1, reaction_number=1)

        self.assertEqual(result, [46.07])
