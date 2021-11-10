
#define interrupt_pin 3 // connect DAQ output
#define analog_pin 0 // connect DAQ output (same as interrupt_pin)

#define stepper_dir_pin 2 // connect to DIR terminal
#define stepper_step_pin 9 // connect to STEP terminal

#define step_delay 100

void setup() {
  // Stepper control pins
  pinMode(stepper_dir_pin, OUTPUT);
  pinMode(stepper_step_pin, OUTPUT);

  // DAQ read interrupt pin
  attachInterrupt(digitalPinToInterrupt(interrupt_pin), step_motor, CHANGE);
}

void loop() {
  // just wait for interrupts
}

// Step the motor in correct direction based on value of analog pin
void step_motor() {
  A = analogRead(analog_pin);  // read the input pin
  V = map(A, 0, 1023, 0, 5); // map to voltage

  if (V > 4) { // step CW
    digitalWrite(stepper_dir_pin, HIGH); // CW
    step(step_delay)
  }
  else if (V < 1) { // step CCW
    digitalWrite(stepper_dir_pin, LOW); // CCW
    step(step_delay)

  else { // no step
        
  }
}

// Step the motor
void step(delay) {  
  digitalWrite(stepper_step_pin, HIGH);
  delayMicroseconds(delay);
  digitalWrite(stepper_step_pin, LOW);
}