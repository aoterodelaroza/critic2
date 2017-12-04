  RegisterCallback(BIND_QUIT,(void *) quit_callback,(void *) rootwin);

static void *eventbind[BIND_MAX]; // bind -> event
static void *databind[BIND_MAX]; // bind -> data

    eventbind[i] = nullptr;
    databind[i] = nullptr;

void RegisterCallback(int bind,void *callback,void *data){
  eventbind[bind] = callback;
  databind[bind] = data;
}

void ProcessCallbacks(){
  ImGuiIO& io = GetIO();

  for (auto it=keymap.begin(); it != keymap.end(); it++){
    if (eventbind[it->second] && IsKeyPressed(it->first.first,false) && IsModPressed(it->first.second)){
      ((void (*)(void *)) eventbind[it->second])(databind[it->second]);
    }
  }
}

