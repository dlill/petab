# def addEvent(self, trigger, assignments, persistent=True, \
#              initial_value=False, priority=0, delay=0, event_id=''):           
#   """Adds an `Event <http://sbml.org/Software/libSBML/docs/python-api/classlibsbml_1_1_event.html>`_
#         to the model.
# 
#         *trigger* is the string representation of a
#         logical expression that defines when an event is 'triggered', meaning
#         when the event is ready to be executed.
# 
#         *delay* is a numerical value that defines the amount of time between when
#         the event is triggered and when the event assignment is implemented, in
#         previously defined model-wide time units.
# 
#         *assignments* is a dictionary
#         where the keys are variables to be changed and the values are the
#         variables' new values.
# 
#         *persistent* is a boolean that defines whether the
#         event will still be executed if the trigger switches from ``True`` to ``False``
#         between the event's trigger and its execution.
# 
#         *initial_value* is the value of *trigger* when t < 0.
# 
#         *priority* is a numerical value that determines
#         which event is executed if two events are executed at the same time.  The
#         event with the larger ``priority`` is executed.
# 
#         .. note:: An event is only triggered when the trigger switches from ``False`` to
#             ``True``.  If the trigger's initial value is ``True``, the event will not be
#             triggered until the value switches to ``False`` and then back to ``True``."""
# 
# e1 = self.model.createEvent()
# self._check(e1,     'create event')
# if len(event_id) == 0:
#   event_id = 'e' + str(self.model.getNumEvents())
# self._check(e1.setId(event_id),    'add id to event')
# if self.document.getLevel()==3 or (self.document.getLevel()==2 \
#                                    and self.document.getVersion()==4):
#   self._check(e1.setUseValuesFromTriggerTime(True), 'set use values from trigger time')
# 
# tri = e1.createTrigger()
# self._check(tri,  'add trigger to event')
# tri_ast = libsbml.parseL3Formula(trigger)
# self._check(tri.setMath(tri_ast),     'add formula to trigger')
# if self.document.getLevel() == 3:
#   self._check(tri.setPersistent(persistent),   'set persistence of trigger')
# self._check(tri.setInitialValue(initial_value), 'set initial value of trigger')
# 
# de = e1.createDelay()
# if self.document.getLevel() == 3:
#   k = self.addParameter(event_id+'Delay', delay, self.model.getTimeUnits())
# else:
#   k = self.addParameter(event_id+'Delay', delay, 'time')
# self._check(de,               'add delay to event')
# delay_ast = libsbml.parseL3Formula(k.getId())
# self._check(de.setMath(delay_ast),     'set formula for delay')
# 
# for a in assignments.keys():
#   assign = e1.createEventAssignment()
# self._check(assign,   'add event assignment to event')
# self._check(assign.setVariable(a),  'add variable to event assignment')
# val_ast = libsbml.parseL3Formula(assignments.get(a))
# self._check(assign.setMath(val_ast),    'add value to event assignment')
# 
# if self.document.getLevel() == 3:
#   pri = e1.createPriority()
# pri_ast = libsbml.parseL3Formula(str(priority))
# self._check(pri.setMath(pri_ast), 'add priority to event')
# return e1

library(libSBML)
devtools::load_all("~/Promotion/Promotion/Projects/petab")
#' .. Eqnlist and objects -----
equationList <- NULL
equationList <- addReaction(equationList, from = "E + S", to = "ES", rate = "(kon)*E*S",
                  description = "production of complex")
equationList <- addReaction(equationList, from = "ES", to = "E + S", rate = "koff*ES",
                  description = "decay of complex")
equationList <- addReaction(equationList, from = "ES", to = "E + P", rate = "kcat*ES",
                  description = "production of product")
equationList <- eqnlist_addDefaultCompartment(equationList, "cytoplasm")

filename <- modelname <- file.path(tempdir(), "model.xml")

parInfo <- data.table(tibble::tribble(
  ~parName, ~parValue, ~parUnit,
  "kon"   ,     0.10,"litre_per_mole_per_second" ,
  "koff"  ,     0.55,"per_second" ,
  "kcat"  ,     1.00,"per_second" ))



unitInfo <- getUnitInfo()
compartmentInfo <- getCompartmentInfo(equationList)
speciesInfo <- getSpeciesInfo(equationList)

eventList <- eventlist(var = "E", time = 0, value = 1, root = NA, method = "replace")
eventList <- addEvent(eventList, var = "E", time = 1, value = 0, root = NA, method = "replace")
eventInfo <- getEventInfo(eventList)


# ..  -----

# Collect arguments
reactionInfo <- getReactionInfo(equationList,parInfo = parInfo)

unitInfoList        <- purrr::transpose(unitInfo)
speciesInfoList     <- purrr::transpose(speciesInfo)
parInfoList         <- purrr::transpose(parInfo)
compartmentInfoList <- purrr::transpose(compartmentInfo)
reactionInfoList    <- purrr::transpose(reactionInfo)
eventInfoList    <- purrr::transpose(eventInfo)

# Start SBML document
sbmlDoc = SBMLDocument(level = 3, version = 2) 
model <- sbml_initialize(modelname, sbmlDoc)

# Populate with content
for (x in unitInfoList)        do.call(sbml_addOneUnit,        c(list(model = model),x))
for (x in compartmentInfoList) do.call(sbml_addOneCompartment, c(list(model = model),x))
for (x in speciesInfoList)     do.call(sbml_addOneSpecies,     c(list(model = model),x))
for (x in parInfoList)         do.call(sbml_addOneParameter,   c(list(model = model),x))
for (x in reactionInfoList)    do.call(sbml_addOneReaction,    c(list(model = model),x))
for (x in eventInfoList)       do.call(sbml_addOneEvent,       c(list(model = model),x))




# ..  -----
