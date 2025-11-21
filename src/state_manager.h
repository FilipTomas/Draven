
#ifndef RAVEN_STATE_MANAGER_H_
#define RAVEN_STATE_MANAGER_H_

namespace raven {

enum struct GraphState {
  Construct_Graph,
  Construct_Graph_2,
  Assemble_Transitive_Edges,
  Assemble_Tips_Bubbles,
  Assemble_Force_Directed,
  Assemble_Diploids
};

class StateManager {
public:
  StateManager() = default;

  StateManager(GraphState initialState) : state_(initialState) {}

  GraphState state() const {
    return state_;
  }

  void set_state(GraphState newState) {
    state_ = newState;
  }

  void advance_state() {
    state_ = (GraphState)((int)state_ + 1);
  }

  bool construct_any_state() {
    if (state_ == GraphState::Construct_Graph
      || state_ == GraphState::Construct_Graph_2) return true;

    return false;
  }

  bool assemble_any_state() {
    if (state_ == GraphState::Assemble_Transitive_Edges
      || state_ == GraphState::Assemble_Tips_Bubbles
      || state_ == GraphState::Assemble_Force_Directed) return true;

    return false;
  }

private:
  GraphState state_;
};

}

#endif //RAVEN_STATE_MANAGER_H_
