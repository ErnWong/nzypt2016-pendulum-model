from kivy.app import App
from kivy.uix.widget import Widget
from kivy.properties import NumericProperty
from kivy.clock import Clock
import sys

class PendulumScreen(Widget):
    simtime = NumericProperty(100.0)
    pivotx = NumericProperty(100.0)
    pivoty = NumericProperty(100.0)
    pivotz = NumericProperty(100.0)
    bobx = NumericProperty(100.0)
    boby = NumericProperty(100.0)
    bobz = NumericProperty(100.0)
    def update(self, dt):
        line = sys.stdin.readline()
        self.simtime, self.pivotx,self.pivoty,self.pivotz,self.bobx,self.boby,self.bobz = map(float, line.split())
        print('time:' + str(self.simtime) + ' | x:' + str(self.bobx) + ' | y:' + str(self.boby) + ' | z:' + str(self.bobz))

class PendulumApp(App):
    def build(self):
        screen = PendulumScreen()
        Clock.schedule_interval(screen.update, 1.0/60.0)
        return screen

if __name__ == '__main__':
    PendulumApp().run()
