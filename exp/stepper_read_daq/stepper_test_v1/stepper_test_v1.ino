/* Example sketch to control a stepper motor with TB6600 stepper motor driver and Arduino without a library: number of revolutions, speed and direction. More info: https://www.makerguides.com */

// Define stepper motor connections and steps per revolution:
#define dirPin 10
#define stepPin 9
#define step_size 1.8
#define micro_step 16
#define step_delay 100
#define period_delay 50

int A = 100;
float res = step_size / micro_step;
int stepA = round(A /  res);

void setup() {
  // Declare pins as output:
  pinMode(stepPin, OUTPUT);
  pinMode(dirPin, OUTPUT);
}

void loop() {
  // Set the spinning direction clockwise:
  digitalWrite(dirPin, HIGH);

  // Spin the stepper motor 1 revolution slowly:
  for (int i = 0; i < stepA; i++) {
    // These four lines result in 1 step:
    digitalWrite(stepPin, HIGH);
    delayMicroseconds(step_delay);
    digitalWrite(stepPin, LOW);
    //delayMicroseconds(step_delay);
  }
    
  delay(period_delay/2);

  // Set the spinning direction clockwise:
  digitalWrite(dirPin, LOW);

  // Spin the stepper motor 1 revolution slowly:
  for (int i = 0; i < stepA; i++) {
    // These four lines result in 1 step:
    digitalWrite(stepPin, HIGH);
    delayMicroseconds(step_delay);
    digitalWrite(stepPin, LOW);
    //delayMicroseconds(step_delay);
  }    
    
  delay(period_delay/2);

}