"""Explicit route validation with compatibility for existing modality/use_rfd3 inputs."""

SMALL_MOLECULE = "small_molecule"
PROTEIN_BINDER = "protein_binder_rfd3"
ADC = "adc"
DEGRADER = "degrader"
MODALITY_ROUTES = {"SMALL_MOLECULE": SMALL_MOLECULE, "DE_NOVO_BINDER": PROTEIN_BINDER,
                   "ADC": ADC, "DEGRADER": DEGRADER}
ROUTES = set(MODALITY_ROUTES.values())


def validate_route(route):
    route = MODALITY_ROUTES.get(route, PROTEIN_BINDER if route == "de_novo_binder" else route)
    if route not in ROUTES:
        raise ValueError(f"Unknown route: {route!r}")
    return route


def normalize_route_state(state):
    state = dict(state)
    if "route" not in state:
        return state  # Legacy graph/state APIs remain supported.
    route = validate_route(state["route"])
    state["route"] = route
    if route in {ADC, DEGRADER}:
        state.update(design_mode=route, use_rfd3=False)
        return state
    expected = route == PROTEIN_BINDER
    if "use_rfd3" in state and state["use_rfd3"] is not expected:
        raise ValueError("route conflicts with use_rfd3")
    state.update(
        use_rfd3=expected,
        design_mode="protein_binder" if expected else "small_molecule",
    )
    return state


def validate_handoff_routes(designs, requested=None):
    """Check the complete input before either branch can start a backend."""
    for design in designs:
        modality = design.get("modality")
        expected = MODALITY_ROUTES.get(modality)
        if expected is None:
            raise ValueError(f"Unsupported modality: {modality!r}")
        route = design.get("route", expected)
        if route is not None:
            validate_route(route)
            if expected and route != expected:
                raise ValueError("Design route conflicts with modality")
            if requested and route != requested:
                raise ValueError("Design route conflicts with requested route")
        for candidate in design.get("candidates", design.get("ligands", [])):
            if "route" in candidate and validate_route(candidate["route"]) != route:
                raise ValueError("Candidate route conflicts with Design route")
