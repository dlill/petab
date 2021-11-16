each rds contains one or multiple parameterSets

* parameterSet: one parameter vector, i.e. one entry
* parameterSetId
    * Way to navigate to the parameterSet in an unambiguous manner
    * e.g. object_subclass1_subclass2_...
        * mstrust_fitrank1 (might be ambiguous if one adds more fits later)
        * or mstrust_step1_first
        * profileParameterSetId_parameter_optimum
        * profileParameterSetId_parameter_leftEndpoint
        * generic: parframeId_row


* Need to agree on names
* Proposal
    * One unambiguous way, immutable (does not change if fits are added)
    * Implement shortcuts
